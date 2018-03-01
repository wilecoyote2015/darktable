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
  float neighborhood_size;
  float patch_size;
  float h;
} dt_iop_rawdenoise_nlmeans_params_t;

typedef struct dt_iop_rawdenoise_nlmeans_gui_data_t
{
  GtkWidget *stack;
  GtkWidget *box_raw;
  GtkWidget *neighborhood_size;
  GtkWidget *patch_size;
  GtkWidget *h;
  GtkWidget *label_non_raw;
} dt_iop_rawdenoise_nlmeans_gui_data_t;

typedef struct dt_iop_rawdenoise_nlmeans_data_t
{
  float neighborhood_size;
  float patch_size;
  float h;
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
  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "filter strength"));
//  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "patch size"));
//  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "neighborhood size"));
}

void connect_key_accels(dt_iop_module_t *self)
{
  dt_iop_rawdenoise_nlmeans_gui_data_t *g = (dt_iop_rawdenoise_nlmeans_gui_data_t *)self->gui_data;

  dt_accel_connect_slider_iop(self, "filter_strength", GTK_WIDGET(g->h));
  dt_accel_connect_slider_iop(self, "patch size", GTK_WIDGET(g->patch_size));
  dt_accel_connect_slider_iop(self, "neighborhood size", GTK_WIDGET(g->neighborhood_size));
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

static inline void transform_bayer(float *const buf, const int width, const int height, const float a[4],
                                       const float b[4])
{
  // transform pixels for each color sepatarely
  for (int y = 0; y < height; y += 2)
  {
    for (int x = 0; x < width; x+= 2)
    {
      // todo: take care of orders in elements in a and b according to bayer pattern!
      for (int color_y = 0; color_y < 2; color_y++)
      {
        for (int color_x = 0; color_x < 2; color_x++)
        {
          int index_image = (y + color_y) * width + x + color_x;
          int index_color = 2 * color_y + color_x;

          buf[index_image] = 2. * sqrtf((buf[index_image] - b[index_color])/ a[index_color] + 3. / 8.);
        }
      }
    }
  }
}

static inline void backtransform_bayer(float *const buf, const int width, const int height, const float a[4],
                                 const float b[4])
{
  // transform pixels for each color sepatarely
  for (int y = 0; y < height; y += 2)
  {
    for (int x = 0; x < width; x+= 2)
    {
      // todo: take care of orders in elements in a and b according to bayer pattern!
      for (int color_y = 0; color_y < 2; color_y++)
      {
        for (int color_x = 0; color_x < 2; color_x++)
        {
          int index_image = (y + color_y) * width + x + color_x;
          int index_color = 2 * color_y + color_x;

          float value = buf[index_image];
          if (value > 0) {
            float value_tranformed = 1. / 4. * value + 1. / 4. * sqrtf(3. / 2.) / value - 11. / 8. * 1.0 / (value * value)
                                     + 5. / 8. * sqrtf(3. / 2.) * 1.0 / (value * value * value) - 1. / 8.;
            buf[index_image] = a[index_color] * value_tranformed + b[index_color];
          } else {
            buf[index_image] = b[index_color];
          }
        }
      }
    }
  }
}

static inline int index_coords(const int x, const int y, const int width)
{
  return y * width + x;  // todo: correct?
}

void apply_nlmeans(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
                   void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{

  // todo: the input and output are now uint16! Modify the code accordingly!
  // todo: handle divisions by zero, also in anscombe transforms!

  // this is called for preview and full pipe separately, each with its own pixelpipe piece.
  // get our data struct:
  const struct dt_iop_rawdenoise_nlmeans_params_t *const d = (const dt_iop_rawdenoise_nlmeans_params_t *const)piece->data;

  // TODO: fixed K to use adaptive size trading variance and bias!
  // adjust to zoom size:
  const int patch_size = (int) d->patch_size; // pixel filter size
  const int neighborhood_size = (int)d->neighborhood_size;         // nbhood
  const float h = d->h;
  const int num_pixels_patch = (patch_size * 2 + 1) * (patch_size * 2 + 1);

  // todo: there should be given by profile
  const float aa[4] = {0.9f, 0.9f, 0.9f, 0.9f};
  const float bb[4] = {-150.0f, -150.0f, -150.0f, -150.0f};

  const int size_raw_pattern = 2; // todo: derive from sensor type

  const int width = roi_in->width;
  const int height = roi_in->height;

  // get maximum and minimum in dataset // todo: only for testing!
  float minimum = 100.0f, maximum = 0.0f;
  float* p_in = (float*) ivoid;
  for (int index = 0; index < width * height; index++)
  {
    float value = p_in[index];
    if (value < minimum) minimum = value;
    if (value > maximum) maximum = value;
  }



  // P == 0 : this will degenerate to a (fast) bilateral filter.

  float *square_differences = dt_alloc_align(64, (size_t) sizeof(float) * width * height);
  float *weigths_summed = dt_alloc_align(64, (size_t) sizeof(float) * width * height);
  // we want to sum up weights in col[1], and sum up results in col[0], so need to init to 0:
  // output
  memset(ovoid, 0x0, (size_t)sizeof(float) * width * height);

  // input. will be filled by precondition by anscombe transformed data
  float *in = dt_alloc_align(64, (size_t) sizeof(float) * width * height);

  // transform data // todo still wrong
//  transform_bayer(in, width, height, aa, bb);


  // because we don't perdfrm anscombe transform yet, simply copy input to in
  memcpy(in, ivoid, (size_t)sizeof(float) * width * height);

  // get float pointer from ovoid
  float *const out = ((float *const)ovoid);

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
  for(int shift_y = -neighborhood_size_scaled; shift_y <= neighborhood_size_scaled; shift_y += size_raw_pattern)
  {
    for(int shift_x = -neighborhood_size_scaled; shift_x <= neighborhood_size_scaled; shift_x += size_raw_pattern)
    {
      // calculate square differences for current shift
      // todo: in all the for-loops, swap x and y loops maybe, and pre-calculate x*width. We don't want to calc it every time!
      // todo: in some cases, we don't need to calc coordinates, but do that in the for loops!
      for(int y = 0, y_shifted = shift_y; y < height; y++, y_shifted++)
      {
        if (y_shifted < 0 || y_shifted >= height) continue;

        int index_at_y = index_coords(0, y, width);
        int index_at_y_shifted = index_coords(0, y_shifted, width);

        for(int x = 0, x_shifted = shift_x; x < width; x++, x_shifted++)
        {
          if (x_shifted < 0 || x_shifted >= width) continue;
          int index_at_xy = index_at_y + x;
          float difference = in[index_at_xy] - in[index_at_y_shifted + x_shifted];
          square_differences[index_at_xy] = difference * difference;
        }
      }

      // for each pixel, calculate the summed quare difference it's patch and the weight from that.
      // add the value of shifted center pixel, multiplied with weight, to the output, and add weight to summed weights
      for(int y = 0, y_shifted = shift_y; y < height; y++, y_shifted++)
      {
        if (y_shifted < 0 || y_shifted >= height) continue;

        int index_at_y = index_coords(0, y, width);
        int index_at_y_shifted = index_coords(0, y_shifted, width);

        for(int x = 0, x_shifted = shift_x; x < width; x++, x_shifted++)
        {
          if (x_shifted < 0 || x_shifted >= width) continue;
          // shifted indices are shifted patch centers

          // todo: move condition calcs out of loop
          // calculate distance for patch
          float distance_patch = 0.0f;
          for(int y_patch = y - patch_size; y_patch <= y + patch_size; y_patch++)
          {
            if (y_patch < 0 || y_patch >= height) continue;
            int index_at_y_patch = index_coords(0, y_patch, width);
            for(int x_patch = x - patch_size; x_patch <= x + patch_size; x_patch++)
            {
              if (x_patch < 0 || x_patch >= width) continue;
              distance_patch += square_differences[index_at_y_patch + x_patch];
            }
          }

          // normalize distance by patch size
          distance_patch /= (float) num_pixels_patch;

          // calculate weight
          float weight = calculate_weight(distance_patch, h);
          int index_at_xy = index_at_y + x;
          weigths_summed[index_at_xy] += weight;
          out[index_at_xy] += in[index_at_y_shifted + x_shifted] * weight;

        }
      }
    }
  }

  // normalize
  // todo: do coordinate calcs only one time, or better, make it in the for lop head!
  for(int y = 0; y < height; y++) {
    int index_at_y = index_coords(0, y, width);
    for (int x = 0; x < width; x++)
    {
      int index_at_xy = index_at_y + x;
      out[index_at_xy] /= weigths_summed[index_at_xy];
    }
  }


  // free shared tmp memory:
  // todo: re-enable, but here is a crash!
//  dt_free_align(in);
//  dt_free_align(square_differences);
//  dt_free_align(weigths_summed);

  // reverse anscombe transform and scaling to data space
  // todo!
//  backtransform_bayer((float *)ovoid, width, height, aa, bb);

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
  dt_iop_rawdenoise_nlmeans_params_t tmp = (dt_iop_rawdenoise_nlmeans_params_t){ .h = 1.0f , .neighborhood_size = 4.0f,
                                                                                  .patch_size = 4.0f};
//  dt_iop_rawdenoise_nlmeans_params_t tmp = (dt_iop_rawdenoise_nlmeans_params_t){ .neighborhood_size = 4 };
//  dt_iop_rawdenoise_nlmeans_params_t tmp = (dt_iop_rawdenoise_nlmeans_params_t){ .patch_size = 4 };

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
  module->priority = 13; // module order created by iop_dependencies.py, do not edit!
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

  d->patch_size = (int)p->patch_size;
  d->neighborhood_size = (int)p->neighborhood_size;
  d->h = p->h;


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

  dt_bauhaus_slider_set(g->patch_size, p->patch_size);
  dt_bauhaus_slider_set(g->neighborhood_size, p->neighborhood_size);
  dt_bauhaus_slider_set(g->h, p->h);

  gtk_stack_set_visible_child_name(GTK_STACK(g->stack), self->hide_enable_button ? "non_raw" : "raw");
}

static void patch_size_callback(GtkWidget *slider, gpointer user_data)
{
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  if(self->dt->gui->reset) return;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;
  p->patch_size = (int)dt_bauhaus_slider_get(slider);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void neighborhood_size_callback(GtkWidget *slider, gpointer user_data)
{
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  if(self->dt->gui->reset) return;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;
  p->neighborhood_size = (int)dt_bauhaus_slider_get(slider);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void h_callback(GtkWidget *slider, gpointer user_data)
{
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  if(self->dt->gui->reset) return;
  dt_iop_rawdenoise_nlmeans_params_t *p = (dt_iop_rawdenoise_nlmeans_params_t *)self->params;
  p->h = dt_bauhaus_slider_get(slider);
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

  // neighborhood_size
  g->neighborhood_size = dt_bauhaus_slider_new_with_range(self, 1.0f, 10.0f, 1., p->neighborhood_size, 0);
  gtk_box_pack_start(GTK_BOX(g->box_raw), GTK_WIDGET(g->neighborhood_size), TRUE, TRUE, 0);
  dt_bauhaus_widget_set_label(g->neighborhood_size, NULL, _("neighborhood size"));
  g_signal_connect(G_OBJECT(g->neighborhood_size), "value-changed", G_CALLBACK(neighborhood_size_callback), self);

  // patch_size
  g->patch_size = dt_bauhaus_slider_new_with_range(self, 1.0f, 10.0f, 1., p->patch_size, 0);
  gtk_box_pack_start(GTK_BOX(g->box_raw), GTK_WIDGET(g->patch_size), TRUE, TRUE, 0);
  dt_bauhaus_widget_set_label(g->patch_size, NULL, _("patch size"));
  g_signal_connect(G_OBJECT(g->patch_size), "value-changed", G_CALLBACK(patch_size_callback), self);

  // h
  g->h = dt_bauhaus_slider_new_with_range(self, 0.01f, 2.0f, 0.01, p->h, 2);
  gtk_box_pack_start(GTK_BOX(g->box_raw), GTK_WIDGET(g->h), TRUE, TRUE, 0);
  dt_bauhaus_widget_set_label(g->h, NULL, _("filter strength"));
  g_signal_connect(G_OBJECT(g->h), "value-changed", G_CALLBACK(h_callback), self);

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
