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
  const float exponent = value / h*h;
  return fast_mexp2f(exponent);
}

void apply_nlmeans(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
                   void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  // this is called for preview and full pipe separately, each with its own pixelpipe piece.
  // get our data struct:
  const dt_iop_rawdenoise_nlmeans_params_t *const d = (dt_iop_rawdenoise_nlmeans_params_t *)piece->data;



//  // adjust to zoom size:
//  const int P = ceilf(d->radius * fmin(roi_in->scale, 2.0f) / fmax(piece->iscale, 1.0f)); // pixel filter size
//  const int K = ceilf(21 * fmin(roi_in->scale, 2.0f) / fmax(piece->iscale, 1.0f));         // nbhood
//  const float sharpness = 3000.0f / (1.0f + d->strength);
//  if(P < 1)
//  {
//    // nothing to do from this distance:
//    memcpy(ovoid, ivoid, (size_t)sizeof(float) * 4 * roi_out->width * roi_out->height);
//    return;
//  }

  // patch and neighorhood size
  // todo: get from data
  const int size_patch = 7;
  const int size_neighborhood = 10;
  const float h = 1.0f;  // todo: probably wrong value


  // get norm to norm data to 1
  // todo: this will not be necessary with variance stabilization
  float max_value = 1.0f; // todo: this is most probably not correct!
  const float norm = 1.0f / max_value;

  float *Sa = dt_alloc_align(64, (size_t)sizeof(float) * roi_out->width * dt_get_num_threads());
  // we want to sum up weights in col[3], so need to init to 0:
  memset(ovoid, 0x0, (size_t)sizeof(float) * roi_out->width * roi_out->height);  // todo: originally was *4, but we dont have multiple colors, so i guess it must be omitted

  // for each shift vector
  // todo: adjust shift and max value to raw pattern size
  for(int kj = -size_neighborhood; kj <= size_neighborhood; kj++)
  {
    for(int ki = -size_neighborhood; ki <= size_neighborhood; ki++)
    {
      int inited_slide = 0;
// don't construct summed area tables but use sliding window! (applies to cpu version res < 1k only, or else
// we will add up errors)
// do this in parallel with a little threading overhead. could parallelize the outer loops with a bit more
// memory
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) firstprivate(inited_slide) shared(kj, ki, Sa)
#endif
      for(int j = 0; j < roi_out->height; j++)
      {
        // if shifted position is out of image range, skip
        if(j + kj < 0 || j + kj >= roi_out->height) continue;

        float *S = Sa + (size_t)dt_get_thread_num() * roi_out->width;

        // pointer to input and output pixel column
        const float *ins = ((float *)ivoid) + ((size_t)roi_in->width * (j + kj) + ki);
        float *out = ((float *)ovoid) + (size_t)roi_out->width * j;

        const int Pm = MIN(MIN(size_patch, j + kj), j);
        const int PM = MIN(MIN(size_patch, roi_out->height - 1 - j - kj), roi_out->height - 1 - j);
        // first line of every thread
        // TODO: also every once in a while to assert numerical precision!
        if(!inited_slide)
        {
          // sum up a line
          memset(S, 0x0, sizeof(float) * roi_out->width);
          for(int jj = -Pm; jj <= PM; jj++)
          {
            int i = MAX(0, -ki);
            float *s = S + i;
            const float *inp = ((float *)ivoid) + i + (size_t)roi_in->width * (j + jj);
            const float *inps = ((float *)ivoid) + i + ((size_t)roi_in->width * (j + jj + kj) + ki);
            const int last = roi_out->width + MIN(0, -ki);
            for(; i < last; i++, inp ++, inps++, s++)
            {
              // todo: inp and inps were array elements with indices over color channels. now are simpy value references
              *s += (*inp - *inps) * (*inp - *inps) * norm;
            }
          }
          // only reuse this if we had a full stripe
          if(Pm == size_patch && PM == size_patch) inited_slide = 1;
        }

        // sliding window for this line:
        float *s = S;
        float slide = 0.0f;
        // sum up the first -P..P
        for(int i = 0; i < 2 * size_patch + 1; i++) slide += s[i];
        for(int i = 0; i < roi_out->width; i++, s++, ins++, out++)
        {
          if(i - size_patch > 0 && i + size_patch < roi_out->width) slide += s[size_patch] - s[-size_patch - 1];
          if(i + ki >= 0 && i + ki < roi_out->width)
          {
            // todo: what is iv?
            const float iv = *ins;

            // calculate the value of nl mean using weight
            *out += iv * calculate_weight(slide, h);
          }
        }
        if(inited_slide && j + size_patch + 1 + MAX(0, kj) < roi_out->height)
        {
          // sliding window in j direction:
          int i = MAX(0, -ki);
          s = S + i;
          const float *inp = ((float *)ivoid) + i + (size_t)roi_in->width * (j + size_patch + 1);
          const float *inps = ((float *)ivoid) + i + ((size_t)roi_in->width * (j + size_patch + 1 + kj) + ki);
          const float *inm = ((float *)ivoid) + i + (size_t)roi_in->width * (j - size_patch);
          const float *inms = ((float *)ivoid) + i + ((size_t)roi_in->width * (j - size_patch + kj) + ki);
          const int last = roi_out->width + MIN(0, -ki);
          for(; i < last; i++, inp++, inps++, inm++, inms++, s++)
          {
            float stmp = s[0];
            stmp += ((*inp - *inps) * (*inp - *inps) - (*inm - *inms) * (*inm - *inms)) * norm;
            s[0] = stmp;
          }
        }
        else
          inited_slide = 0;
      }
    }
  }

  // todo: normalization for weights?

  // free shared tmp memory:
  dt_free_align(Sa);

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
