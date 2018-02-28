/*
    This file is part of darktable,
    copyright (c) 2011 bruce guenter
    copyright (c) 2012 henrik andersson


    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "control/control.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "iop/iop_api.h"

#include <gtk/gtk.h>
#include <stdlib.h>
#include <strings.h>

DT_MODULE_INTROSPECTION(1, dt_iop_rawdenoise_nlmeans_params_t)

typedef struct dt_iop_rawdenoise_nlmeans_params_t
{
  float threshold;
} dt_iop_rawdenoise_nlmeans_params_t;

typedef struct dt_iop_rawdenoise_nlmeans_gui_data_t
{
  GtkWidget *stack;
  GtkWidget *box_raw;
  GtkWidget *threshold;
  GtkWidget *label_non_raw;
} dt_iop_rawdenoise_nlmeans_gui_data_t;

typedef struct dt_iop_rawdenoise_nlmeans_data_t
{
  float threshold;
} dt_iop_rawdenoise_nlmeans_data_t;

typedef struct dt_iop_rawdenoise_nlmeans_global_data_t
{
} dt_iop_rawdenoise_nlmeans_global_data_t;

const char *name()
{
  return _("raw denoise nl means");
}

int flags()
{
  return IOP_FLAGS_SUPPORTS_BLENDING;
}

int groups()
{
  return IOP_GROUP_CORRECT;
}

void init_key_accels(dt_iop_module_so_t *self)
{
  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "noise threshold"));
}

void connect_key_accels(dt_iop_module_t *self)
{
  dt_iop_rawdenoise_nlmeans_gui_data_t *g = (dt_iop_rawdenoise_nlmeans_gui_data_t *)self->gui_data;

  dt_accel_connect_slider_iop(self, "noise threshold", GTK_WIDGET(g->threshold));
}

typedef union floatint_t
{
  float f;
  uint32_t i;
} floatint_t;


#define BIT16 65536.0


static inline float fast_mexp2f(const float x)
{
  const float i1 = (float)0x3f800000u; // 2^0
  const float i2 = (float)0x3f000000u; // 2^-1
  const float k0 = i1 + x * (i2 - i1);
  floatint_t k;
  k.i = k0 >= (float)0x800000u ? k0 : 0;
  return k.f;
}

static float calculate_weight(const float value, const float h)
{
  const float exponent = value / (h*h);
  return fast_mexp2f(exponent);
}

void apply_nlmeans(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
                   void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  // this is called for preview and full pipe separately, each with its own pixelpipe piece.
  // get our data struct:
  const struct dt_iop_rawdenoise_nlmeans_params_t *const d = (const dt_iop_rawdenoise_nlmeans_params_t *const)piece->data;

  const int n_colors = piece->colors;  // todo: isn't it simply 1?

  // TODO: fixed K to use adaptive size trading variance and bias!
  // adjust to zoom size:
  const float scale = fmin(roi_in->scale, 2.0f) / fmax(piece->iscale, 1.0f);
  const int patch_size = 1; // pixel filter size
  const int neighborhood_size = 1;         // nbhood

  const int size_raw_pattern = 2; // todo: derive from sensor type


  const int n_channels = 2;  // number of channels is number of colors + 1 for weights

  // P == 0 : this will degenerate to a (fast) bilateral filter.

  float *values_summed_all = dt_alloc_align(64, (size_t)sizeof(float) * roi_out->width * dt_get_num_threads());
  // we want to sum up weights in col[1], and sum up results in col[0], so need to init to 0:
  // output
  memset(ovoid, 0x0, (size_t)sizeof(float) * roi_out->width * roi_out->height * n_channels);

  // input. will be filled by precondition by anscombe transformed data
  float *in = dt_alloc_align(64, (size_t) n_channels * sizeof(float) * roi_in->width * roi_in->height);


  // because we don't perdfrm anscombe transform yet, simply copy input to in
  memcpy(in, ivoid, (size_t)sizeof(float) * n_channels * roi_out->width * roi_out->height);

  // transform input data to gaussian noise with std. dev 1
  // todo!
//  const float wb[3] = { piece->pipe->dsc.processed_maximum[0] * d->strength * (scale * scale),
//                        piece->pipe->dsc.processed_maximum[1] * d->strength * (scale * scale),
//                        piece->pipe->dsc.processed_maximum[2] * d->strength * (scale * scale) };
//  const float aa[3] = { d->a[1] * wb[0], d->a[1] * wb[1], d->a[1] * wb[2] };
//  const float bb[3] = { d->b[1] * wb[0], d->b[1] * wb[1], d->b[1] * wb[2] };
//  precondition((float *)ivoid, in, roi_in->width, roi_in->height, aa, bb);

  // for each shift vector

  const int neighborhood_size_scaled = neighborhood_size * size_raw_pattern;
  for(int kj = -neighborhood_size_scaled; kj <= neighborhood_size_scaled; kj += size_raw_pattern)
  {
    for(int ki = -neighborhood_size_scaled; ki <= neighborhood_size_scaled; ki += size_raw_pattern)
    {
      // TODO: adaptive K tests here!
      // TODO: expf eval for real bilateral experience :)

      int inited_slide = 0;
// don't construct summed area tables but use sliding window! (applies to cpu version res < 1k only, or else
// we will add up errors)
// do this in parallel with a little threading overhead. could parallelize the outer loops with a bit more
// memory
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) firstprivate(inited_slide) shared(kj, ki, in, values_summed_all)
#endif
      for(int j = 0; j < roi_out->height; j++)
      {
        if(j + kj < 0 || j + kj >= roi_out->height) continue;
        float *values_summed_line = values_summed_all + dt_get_thread_num() * roi_out->width;
        const float *input_at_j = in + n_channels * ((size_t)roi_in->width * (j + kj) + ki);
        float *out = ((float *)ovoid) + (size_t) n_channels * roi_out->width * j;

        const int Pm = MIN(MIN(patch_size, j + kj), j);
        const int PM = MIN(MIN(patch_size, roi_out->height - 1 - j - kj), roi_out->height - 1 - j);
        // first line of every thread
        // TODO: also every once in a while to assert numerical precision!
        if(!inited_slide)
        {
          // sum up a line
          memset(values_summed_line, 0x0, sizeof(float) * roi_out->width);
          for(int jj = -Pm; jj <= PM; jj++)
          {
            int i = MAX(0, -ki);
            float *s = values_summed_line + i;
            const float *inp = in + n_channels * i + (size_t) n_channels * roi_in->width * (j + jj);
            const float *inps = in + n_channels * i + n_channels * ((size_t)roi_in->width * (j + jj + kj) + ki);
            const int last = roi_out->width + MIN(0, -ki);
            for(; i < last; i++, inp += n_channels, inps += n_channels, s++)
            {
              for(int k = 0; k < 1; k++) s[0] += (inp[k] - inps[k]) * (inp[k] - inps[k]);  // todo: replace 1 by n_colors? even need loop here?
            }
          }
          // only reuse this if we had a full stripe
          if(Pm == patch_size && PM == patch_size) inited_slide = 1;
        }

        // sliding window for this line:
        float *values_summed_window = values_summed_line;
        float slide = 0.0f;  // sum of patch
        // sum up the first -P..P
        for(int i = 0; i < 2 * patch_size + 1; i++) slide += values_summed_window[i];
        for(int i = 0; i < roi_out->width; i++, values_summed_window++, input_at_j += n_channels, out += n_channels)
        {
          // FIXME: the comment above is actually relevant even for 1000 px width already.
          // XXX    numerical precision will not forgive us:
          if(i - patch_size > 0 && i + patch_size < roi_out->width)
            slide += values_summed_window[patch_size] - values_summed_window[-patch_size - 1];
          if(i + ki >= 0 && i + ki < roi_out->width)
          {
            // TODO: could put that outside the loop.
            // DEBUG XXX bring back to computable range:
            // norm the distance according to patch size
            // todo: shouldn't it be pixel count, which is square 2 * patch_size + 1? Is done now. see of this works.
            const int pixel_count_patch = (2 * patch_size + 1) * (2 * patch_size + 1);
            const float norm = 1.0f / (float)pixel_count_patch;  // todo: why the 0.15f ? replaced by 1 now.
            const float iv[2] = { input_at_j[0], 1.0f };  // todo: get this dynamically with n_colors?
#if defined(_OPENMP) && defined(OPENMP_SIMD_)
#pragma omp SIMD()
#endif
            // calculate weight and insert weighted sum in valuue iv[0] and weight in iv[1]
            for(size_t c = 0; c < 2; c++)  // todo: use n_channels?
            {
              out[c] += iv[c] * fast_mexp2f(fmaxf(0.0f, slide * norm - 2.0f) / 10000000000000000000.0f);  // todo: try without substraction of 2sigma.
              // todo: don't divide by large number, but by h**2
            }
          }
        }
        if(inited_slide && j + patch_size + 1 + MAX(0, kj) < roi_out->height)
        {
          // sliding window in j direction:
          int i = MAX(0, -ki);
          values_summed_window = values_summed_line + i;
          const float *inp = in + n_channels * i + n_channels * (size_t)roi_in->width * (j + patch_size + 1);
          const float *inps = in + n_channels * i + n_channels * ((size_t)roi_in->width * (j + patch_size + 1 + kj) + ki);
          const float *inm = in + n_channels * i + n_channels * (size_t)roi_in->width * (j - patch_size);
          const float *inms = in + n_channels * i + n_channels * ((size_t)roi_in->width * (j - patch_size + kj) + ki);
          const int last = roi_out->width + MIN(0, -ki);
          for(; i < last; i++, inp += n_channels, inps += n_channels, inm += n_channels, inms += n_channels, values_summed_window++)
          {
            float stmp = values_summed_window[0];
            for(int k = 0; k < 1; k++)  // todo: replace 1 by n_colors? even need loop here?
              stmp += ((inp[k] - inps[k]) * (inp[k] - inps[k]) - (inm[k] - inms[k]) * (inm[k] - inms[k]));
            values_summed_window[0] = stmp;
          }
        }
        else
          inited_slide = 0;
      }
    }
  }

  float *const out = ((float *const)ovoid);

// normalize
#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static)
#endif
  for(size_t k = 0; k < (size_t)n_colors * roi_out->width * roi_out->height; k += n_colors)
  {
    if(out[k + 1] <= 0.0f) continue;
    for(size_t c = 0; c < 1; c++)  // todo: c < num_colors?
    {
      out[k + c] *= (1.0f / out[k + 1]);  // todo: replace 1 by n_colors?
    }
  }

  // free shared tmp memory:
  dt_free_align(values_summed_all);
  dt_free_align(in);

  // reverse anscombe transform and scaling to data space
  // todo!
//  backtransform((float *)ovoid, roi_in->width, roi_in->height, aa, bb);

  if(piece->pipe->mask_display & DT_DEV_PIXELPIPE_DISPLAY_MASK) dt_iop_alpha_copy(ivoid, ovoid, roi_out->width, roi_out->height);
}

void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
             void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_rawdenoise_nlmeans_data_t *d = (dt_iop_rawdenoise_nlmeans_data_t *)piece->data;

  apply_nlmeans(self, piece, ivoid, ovoid, roi_in, roi_out);
}


void reload_defaults(dt_iop_module_t *module)
{
  // init defaults:
  dt_iop_rawdenoise_nlmeans_params_t tmp = (dt_iop_rawdenoise_nlmeans_params_t){ .threshold = 0.01 };

  // we might be called from presets update infrastructure => there is no image
  if(!module->dev) goto end;

  // can't be switched on for non-raw images:
  if(dt_image_is_raw(&module->dev->image_storage))
    module->hide_enable_button = 0;
  else
    module->hide_enable_button = 1;
  module->default_enabled = 0;

end:
  memcpy(module->params, &tmp, sizeof(dt_iop_rawdenoise_nlmeans_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_rawdenoise_nlmeans_params_t));
}

void init(dt_iop_module_t *module)
{
  module->data = NULL;
  module->params = calloc(1, sizeof(dt_iop_rawdenoise_nlmeans_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_rawdenoise_nlmeans_params_t));
  module->default_enabled = 0;

  // raw denoise must come just before demosaicing.
  module->priority = 103; // module order created by iop_dependencies.py, do not edit!
  module->params_size = sizeof(dt_iop_rawdenoise_nlmeans_params_t);
  module->gui_data = NULL;
}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
  free(module->data);
  module->data = NULL;
}

void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *params, dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)params;
  dt_iop_rawdenoise_nlmeans_data_t *d = (dt_iop_rawdenoise_nlmeans_data_t *)piece->data;

  d->threshold = p->threshold;

  if (!(pipe->image.flags & DT_IMAGE_RAW))
    piece->enabled = 0;
}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_rawdenoise_nlmeans_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
  piece->data = NULL;
}

void gui_update(dt_iop_module_t *self)
{
  dt_iop_rawdenoise_nlmeans_gui_data_t *g = (dt_iop_rawdenoise_nlmeans_gui_data_t *)self->gui_data;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;

  dt_bauhaus_slider_set(g->threshold, p->threshold);

  gtk_stack_set_visible_child_name(GTK_STACK(g->stack), self->hide_enable_button ? "non_raw" : "raw");
}

static void threshold_callback(GtkWidget *slider, gpointer user_data)
{
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  if(self->dt->gui->reset) return;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;
  p->threshold = dt_bauhaus_slider_get(slider);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void gui_init(dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_rawdenoise_nlmeans_gui_data_t));
  dt_iop_rawdenoise_nlmeans_gui_data_t *g = (dt_iop_rawdenoise_nlmeans_gui_data_t *)self->gui_data;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;

  self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, 0));

  g->stack = gtk_stack_new();
  gtk_stack_set_homogeneous(GTK_STACK(g->stack), FALSE);
  gtk_box_pack_start(GTK_BOX(self->widget), g->stack, TRUE, TRUE, 0);

  g->box_raw = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);

  /* threshold */
  g->threshold = dt_bauhaus_slider_new_with_range(self, 0.0, 0.1, 0.001, p->threshold, 3);
  gtk_box_pack_start(GTK_BOX(g->box_raw), GTK_WIDGET(g->threshold), TRUE, TRUE, 0);
  dt_bauhaus_widget_set_label(g->threshold, NULL, _("noise threshold"));
  g_signal_connect(G_OBJECT(g->threshold), "value-changed", G_CALLBACK(threshold_callback), self);

  gtk_widget_show_all(g->box_raw);
  gtk_stack_add_named(GTK_STACK(g->stack), g->box_raw, "raw");

  g->label_non_raw = gtk_label_new(_("raw denoising\nonly works for raw images."));
  gtk_widget_set_halign(g->label_non_raw, GTK_ALIGN_START);

  gtk_widget_show_all(g->label_non_raw);
  gtk_stack_add_named(GTK_STACK(g->stack), g->label_non_raw, "non_raw");

  gtk_stack_set_visible_child_name(GTK_STACK(g->stack), self->hide_enable_button ? "non_raw" : "raw");
}

void gui_cleanup(dt_iop_module_t *self)
{
  free(self->gui_data);
  self->gui_data = NULL;
}
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
