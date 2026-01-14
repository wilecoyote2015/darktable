/*
    This file is part of darktable,
    Copyright (C) 2018-2025 darktable developers.

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

/*** DOCUMENTATION
 *
 * This module performs local contrast enhancement in scene-referred linear RGB space.
 *
 * It works by computing two luminance masks:
 * 1. A pixel-wise luminance (unblurred)
 * 2. A smoothed luminance using edge-aware filters (guided filter or EIGF)
 *
 * The difference between these two masks represents the local detail/contrast.
 * The local contrast is then enhanced by scaling this difference and applying
 * it as an exposure correction to each pixel.
 *
 * The module should be placed early in the pipe (before color profile)
 * as it operates on scene-linear RGB data.
 *
 ***/

 // TODO: blending parameter should be sqrt or log or so because it is very sensitive to small changes 

#include "common/extra_optimizations.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/fast_guided_filter.h"
#include "common/eigf.h"
#include "common/luminance_mask.h"
#include "common/opencl.h"
#include "control/conf.h"
#include "control/control.h"
#include "develop/blend.h"
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
#include "develop/imageop_gui.h"
#include "gui/accelerators.h"
#include "gui/draw.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include "iop/iop_api.h"
#include "common/iop_group.h"

#ifdef _OPENMP
#include <omp.h>
#endif


DT_MODULE_INTROSPECTION(3, dt_iop_local_contrast_rgb_params_t)


#define MIN_FLOAT exp2f(-16.0f)

/**
 * Filter types for detail preservation / smoothing.
 * DT_TONEEQ_NONE is intentionally omitted as it produces no blur,
 * which would result in no local contrast extraction.
 **/
typedef enum dt_iop_local_contrast_rgb_filter_t
{
  DT_LC_AVG_GUIDED = 0, // $DESCRIPTION: "averaged guided filter"
  DT_LC_GUIDED,         // $DESCRIPTION: "guided filter"
  DT_LC_AVG_EIGF,       // $DESCRIPTION: "averaged EIGF"
  DT_LC_EIGF            // $DESCRIPTION: "EIGF"
} dt_iop_local_contrast_rgb_filter_t;


typedef struct dt_iop_local_contrast_rgb_params_t
{
  // Local contrast scaling factor
  float detail_scale;   // $MIN: 0.0 $MAX: 5.0 $DEFAULT: 1.5 $DESCRIPTION: "detail boost"
  float fine_scale;     // $MIN: 0.0 $MAX: 5.0 $DEFAULT: 1.0 $DESCRIPTION: "fine detail"
  float micro_scale;    // $MIN: 0.0 $MAX: 5.0 $DEFAULT: 1.0 $DESCRIPTION: "micro detail"
  float global_scale;   // $MIN: 0.0 $MAX: 5.0 $DEFAULT: 1.0 $DESCRIPTION: "global contrast"

  // Masking parameters
  // Blending is log-encoded because changes in small values are more noticeable
  float blending;       // $MIN: 0.01 $MAX: 100.0 $DEFAULT: 12.0 $DESCRIPTION: "feature scale"
  float feathering;     // $MIN: 0.01 $MAX: 10000.0 $DEFAULT: 5.0 $DESCRIPTION: "edges refinement/feathering"

  dt_iop_local_contrast_rgb_filter_t details; // $DEFAULT: DT_LC_EIGF $DESCRIPTION: "feature extractor"
  dt_iop_luminance_mask_method_t method;      // $DEFAULT: DT_TONEEQ_NORM_2 $DESCRIPTION: "luminance estimator"
  int iterations;       // $MIN: 1 $MAX: 20 $DEFAULT: 1 $DESCRIPTION: "filter diffusion"
} dt_iop_local_contrast_rgb_params_t;


typedef struct dt_iop_local_contrast_rgb_data_t
{
  float detail_scale;
  float fine_scale;
  float micro_scale;
  float global_scale;
  float blending, feathering;
  float scale;
  int radius;
  int iterations;
  dt_iop_luminance_mask_method_t method;
  dt_iop_local_contrast_rgb_filter_t details;
} dt_iop_local_contrast_rgb_data_t;


typedef struct dt_iop_local_contrast_rgb_global_data_t
{
  // Reserved for OpenCL kernels
} dt_iop_local_contrast_rgb_global_data_t;


typedef struct dt_iop_local_contrast_rgb_gui_data_t
{
  // Flags
  gboolean mask_display;

  // Buffer dimensions
  int buf_width;
  int buf_height;
  int pipe_order;

  // Hash for cache invalidation
  dt_hash_t ui_preview_hash;
  dt_hash_t thumb_preview_hash;
  size_t full_preview_buf_width, full_preview_buf_height;
  size_t thumb_preview_buf_width, thumb_preview_buf_height;

  // Cached luminance buffers
  float *thumb_preview_buf_pixel;     // pixel-wise luminance (no blur)
  float *thumb_preview_buf_smoothed;  // smoothed luminance
  float *thumb_preview_buf_smoothed_fine;
  float *thumb_preview_buf_smoothed_micro;
  float *full_preview_buf_pixel;
  float *full_preview_buf_smoothed;
  float *full_preview_buf_smoothed_fine;
  float *full_preview_buf_smoothed_micro;

  // Cache validity
  gboolean luminance_valid;

  // GTK widgets
  GtkWidget *detail_scale, *fine_scale, *micro_scale, *global_scale;
  GtkWidget *blending;
  GtkWidget *method;
  GtkWidget *details, *feathering, *iterations;
  GtkWidget *show_luminance_mask;
} dt_iop_local_contrast_rgb_gui_data_t;


const char *name()
{
  return _("local contrast rgb");
}

const char *aliases()
{
  return _("local contrast|clarity|detail enhancement");
}

const char **description(dt_iop_module_t *self)
{
  return dt_iop_set_description
    (self, _("enhance local contrast by boosting fine details while preserving edges"),
     _("creative"),
     _("linear, RGB, scene-referred"),
     _("linear, RGB"),
     _("linear, RGB, scene-referred"));
}

int default_group()
{
  return IOP_GROUP_BASIC | IOP_GROUP_EFFECTS;
}

int flags()
{
  return IOP_FLAGS_INCLUDE_IN_STYLES | IOP_FLAGS_SUPPORTS_BLENDING;
}

dt_iop_colorspace_type_t default_colorspace(dt_iop_module_t *self,
                                            dt_dev_pixelpipe_t *pipe,
                                            dt_dev_pixelpipe_iop_t *piece)
{
  return IOP_CS_RGB;
}

int legacy_params(dt_iop_module_t *self,
                  const void *const old_params,
                  const int old_version,
                  void **new_params,
                  int32_t *new_params_size,
                  int *new_version)
{
  typedef struct dt_iop_local_contrast_rgb_params_v2_t
  {
    float detail_scale;
    float global_scale;
    float blending;
    float feathering;
    dt_iop_local_contrast_rgb_filter_t details;
    dt_iop_luminance_mask_method_t method;
    int iterations;
  } dt_iop_local_contrast_rgb_params_v2_t;

  typedef struct dt_iop_local_contrast_rgb_params_v1_t
  {
    float detail_scale;
    float blending;
    float feathering;
    dt_iop_local_contrast_rgb_filter_t details;
    dt_iop_luminance_mask_method_t method;
    int iterations;
  } dt_iop_local_contrast_rgb_params_v1_t;

  if(old_version == 1)
  {
    const dt_iop_local_contrast_rgb_params_v1_t *o = (dt_iop_local_contrast_rgb_params_v1_t *)old_params;
    dt_iop_local_contrast_rgb_params_t *n = malloc(sizeof(dt_iop_local_contrast_rgb_params_t));

    n->detail_scale = o->detail_scale;
    n->fine_scale = 1.0f;
    n->micro_scale = 1.0f;
    n->global_scale = 1.0f;
    n->blending = o->blending;
    n->feathering = o->feathering;
    n->details = o->details;
    n->method = o->method;
    n->iterations = o->iterations;

    *new_params = n;
    *new_params_size = sizeof(dt_iop_local_contrast_rgb_params_t);
    *new_version = 3;
    return 0;
  }
  if(old_version == 2)
  {
    const dt_iop_local_contrast_rgb_params_v2_t *o = (dt_iop_local_contrast_rgb_params_v2_t *)old_params;
    dt_iop_local_contrast_rgb_params_t *n = malloc(sizeof(dt_iop_local_contrast_rgb_params_t));

    n->detail_scale = o->detail_scale;
    n->fine_scale = 1.0f;
    n->micro_scale = 1.0f;
    n->global_scale = o->global_scale;
    n->blending = o->blending;
    n->feathering = o->feathering;
    n->details = o->details;
    n->method = o->method;
    n->iterations = o->iterations;

    *new_params = n;
    *new_params_size = sizeof(dt_iop_local_contrast_rgb_params_t);
    *new_version = 3;
    return 0;
  }
  return 1;
}

/**
 * Helper functions
 **/

static void hash_set_get(const dt_hash_t *hash_in,
                         dt_hash_t *hash_out,
                         dt_pthread_mutex_t *lock)
{
  dt_pthread_mutex_lock(lock);
  *hash_out = *hash_in;
  dt_pthread_mutex_unlock(lock);
}


static void invalidate_luminance_cache(dt_iop_module_t *const self)
{
  dt_iop_local_contrast_rgb_gui_data_t *const restrict g = self->gui_data;

  dt_iop_gui_enter_critical_section(self);
  g->luminance_valid = FALSE;
  g->thumb_preview_hash = DT_INVALID_HASH;
  g->ui_preview_hash = DT_INVALID_HASH;
  dt_iop_gui_leave_critical_section(self);
  dt_iop_refresh_all(self);
}


/**
 * Compute pixel-wise luminance mask (no blur)
 **/
__DT_CLONE_TARGETS__
static inline void compute_pixel_luminance_mask(const float *const restrict in,
                                                float *const restrict luminance,
                                                const size_t width,
                                                const size_t height,
                                                const dt_iop_luminance_mask_method_t method)
{
  // No exposure/contrast boost, just compute raw luminance
  luminance_mask(in, luminance, width, height, method, 1.0f, 0.0f, 1.0f);
}


/**
 * Compute smoothed luminance mask using edge-aware filters
 **/
__DT_CLONE_TARGETS__
static inline void compute_smoothed_luminance_mask(const float *const restrict in,
                                                   float *const restrict luminance,
                                                   const size_t width,
                                                   const size_t height,
                                                const dt_iop_local_contrast_rgb_data_t *const d,
                                                const int radius)
{
  // First compute pixel-wise luminance (no boost)
  luminance_mask(in, luminance, width, height, d->method, 1.0f, 0.0f, 1.0f);

  // Then apply the smoothing filter
  switch(d->details)
  {
    case(DT_LC_AVG_GUIDED):
    {
      fast_surface_blur(luminance, width, height, radius, d->feathering, d->iterations,
                        DT_GF_BLENDING_GEOMEAN, d->scale, 0.0f,
                        exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_GUIDED):
    {
      fast_surface_blur(luminance, width, height, radius, d->feathering, d->iterations,
                        DT_GF_BLENDING_LINEAR, d->scale, 0.0f,
                        exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_AVG_EIGF):
    {
      fast_eigf_surface_blur(luminance, width, height,
                             radius, d->feathering, d->iterations,
                             DT_GF_BLENDING_GEOMEAN, d->scale,
                             0.0f, exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_EIGF):
    {
      fast_eigf_surface_blur(luminance, width, height,
                             radius, d->feathering, d->iterations,
                             DT_GF_BLENDING_LINEAR, d->scale,
                             0.0f, exp2f(-14.0f), 4.0f);
      break;
    }
  }
}


/**
 * Apply local contrast enhancement
 *
 * The detail (local contrast) is the log-space difference between pixel luminance
 * and smoothed luminance. Boosting this difference amplifies local details.
 **/
__DT_CLONE_TARGETS__
static inline void apply_local_contrast(const float *const restrict in,
                                        const float *const restrict luminance_pixel,
                                        const float *const restrict luminance_smoothed,
                                        const float *const restrict luminance_smoothed_fine,
                                        const float *const restrict luminance_smoothed_micro,
                                        float *const restrict out,
                                        const dt_iop_roi_t *const roi_in,
                                        const dt_iop_roi_t *const roi_out,
                                        const dt_iop_local_contrast_rgb_data_t *const d)
{
  const size_t npixels = (size_t)roi_in->width * roi_in->height;

  DT_OMP_FOR()
  for(size_t k = 0; k < npixels; k++)
  {
    // Detail in log space (EV): how much brighter/darker is this pixel
    // compared to its local neighborhood
    // detail = log2(pixel_lum / smoothed_lum) = log2(pixel_lum) - log2(smoothed_lum)
    const float lum_pixel = fmaxf(luminance_pixel[k], MIN_FLOAT);
    const float lum_smoothed = fmaxf(luminance_smoothed[k], MIN_FLOAT);
    const float detail_ev = log2f(lum_pixel / lum_smoothed);

    // Scale the detail: detail_scale = 1.0 means no change
    // > 1.0 boosts local contrast, < 1.0 reduces it
    const float scaled_detail_ev = d->detail_scale * detail_ev;

    // The correction is the difference between scaled and original detail
    float correction_ev = scaled_detail_ev - detail_ev;

    if(luminance_smoothed_fine)
    {
      const float lum_smoothed_fine = fmaxf(luminance_smoothed_fine[k], MIN_FLOAT);
      const float detail_ev_fine = log2f(lum_pixel / lum_smoothed_fine);
      const float scaled_detail_ev_fine = d->fine_scale * detail_ev_fine;
      correction_ev += scaled_detail_ev_fine - detail_ev_fine;
    }

    if(luminance_smoothed_micro)
    {
      const float lum_smoothed_micro = fmaxf(luminance_smoothed_micro[k], MIN_FLOAT);
      const float detail_ev_micro = log2f(lum_pixel / lum_smoothed_micro);
      const float scaled_detail_ev_micro = d->micro_scale * detail_ev_micro;
      correction_ev += scaled_detail_ev_micro - detail_ev_micro;
    }

    // Apply correction in linear space
    // global_scale has the same range as detail_scale.
    const float multiplier = exp2f(correction_ev) * powf(lum_smoothed / 0.1845f, d->global_scale) * 0.1845f / lum_smoothed;

    for_each_channel(c)
      out[4 * k + c] = in[4 * k + c] * multiplier;
  }
}


/**
 * Display the detail mask (difference between pixel and smoothed luminance)
 * Output is a grayscale image normalized to [0, 1] where:
 * - 0.5 = no local detail (pixel matches neighborhood)
 * - < 0.5 = pixel darker than neighborhood
 * - > 0.5 = pixel brighter than neighborhood
 **/
__DT_CLONE_TARGETS__
static inline void display_detail_mask(const float *const restrict luminance_pixel,
                                       const float *const restrict luminance_smoothed,
                                       float *const restrict out,
                                       const size_t width,
                                       const size_t height)
{
  const size_t npixels = width * height;

  DT_OMP_FOR()
  for(size_t k = 0; k < npixels; k++)
  {
    const float lum_pixel = fmaxf(luminance_pixel[k], MIN_FLOAT);
    const float lum_smoothed = fmaxf(luminance_smoothed[k], MIN_FLOAT);

    // Detail in log space, mapped to [0, 1] for display
    // Detail range roughly [-2, +2] EV mapped to [0, 1]
    const float detail_ev = log2f(lum_pixel / lum_smoothed);
    const float intensity = fminf(fmaxf(detail_ev / 4.0f + 0.5f, 0.0f), 1.0f);

    // Set all RGB channels to the same intensity (grayscale)
    out[4 * k + 0] = intensity;
    out[4 * k + 1] = intensity;
    out[4 * k + 2] = intensity;
    // Full opacity
    out[4 * k + 3] = 1.0f;
  }
}


/**
 * Main processing function
 **/
__DT_CLONE_TARGETS__
static void local_contrast_process(dt_iop_module_t *self,
                                   dt_dev_pixelpipe_iop_t *piece,
                                   const void *const restrict ivoid,
                                   void *const restrict ovoid,
                                   const dt_iop_roi_t *const roi_in,
                                   const dt_iop_roi_t *const roi_out)
{
  const dt_iop_local_contrast_rgb_data_t *const d = piece->data;
  dt_iop_local_contrast_rgb_gui_data_t *const g = self->gui_data;

  const float *const restrict in = (float *const)ivoid;
  float *const restrict out = (float *const)ovoid;
  float *restrict luminance_pixel = NULL;
  float *restrict luminance_smoothed = NULL;
  float *restrict luminance_smoothed_fine = NULL;
  float *restrict luminance_smoothed_micro = NULL;

  const size_t width = roi_in->width;
  const size_t height = roi_in->height;
  const size_t num_elem = width * height;

  // Get the hash of the upstream pipe to track changes
  const dt_hash_t hash = dt_dev_pixelpipe_piece_hash(piece, roi_out, TRUE);

  // Sanity checks
  if(width < 1 || height < 1) return;
  if(roi_in->width < roi_out->width || roi_in->height < roi_out->height) return;
  if(piece->colors != 4) return;

  // Init the luminance mask buffers
  gboolean cached = FALSE;

  if(self->dev->gui_attached)
  {
    // If the module instance has changed order in the pipe, invalidate caches
    if(g->pipe_order != piece->module->iop_order)
    {
      dt_iop_gui_enter_critical_section(self);
      g->ui_preview_hash = DT_INVALID_HASH;
      g->thumb_preview_hash = DT_INVALID_HASH;
      g->pipe_order = piece->module->iop_order;
      g->luminance_valid = FALSE;
      dt_iop_gui_leave_critical_section(self);
    }

    if(piece->pipe->type & DT_DEV_PIXELPIPE_FULL)
    {
      // Re-allocate buffers if size changed
      if(g->full_preview_buf_width != width || g->full_preview_buf_height != height)
      {
        dt_free_align(g->full_preview_buf_pixel);
        dt_free_align(g->full_preview_buf_smoothed);
        dt_free_align(g->full_preview_buf_smoothed_fine);
        dt_free_align(g->full_preview_buf_smoothed_micro);
        g->full_preview_buf_pixel = dt_alloc_align_float(num_elem);
        g->full_preview_buf_smoothed = dt_alloc_align_float(num_elem);
        g->full_preview_buf_smoothed_fine = dt_alloc_align_float(num_elem);
        g->full_preview_buf_smoothed_micro = dt_alloc_align_float(num_elem);
        g->full_preview_buf_width = width;
        g->full_preview_buf_height = height;
      }

      luminance_pixel = g->full_preview_buf_pixel;
      luminance_smoothed = g->full_preview_buf_smoothed;
      luminance_smoothed_fine = g->full_preview_buf_smoothed_fine;
      luminance_smoothed_micro = g->full_preview_buf_smoothed_micro;
      cached = TRUE;
    }
    else if(piece->pipe->type & DT_DEV_PIXELPIPE_PREVIEW)
    {
      dt_iop_gui_enter_critical_section(self);
      if(g->thumb_preview_buf_width != width || g->thumb_preview_buf_height != height)
      {
        dt_free_align(g->thumb_preview_buf_pixel);
        dt_free_align(g->thumb_preview_buf_smoothed);
        dt_free_align(g->thumb_preview_buf_smoothed_fine);
        dt_free_align(g->thumb_preview_buf_smoothed_micro);
        g->thumb_preview_buf_pixel = dt_alloc_align_float(num_elem);
        g->thumb_preview_buf_smoothed = dt_alloc_align_float(num_elem);
        g->thumb_preview_buf_smoothed_fine = dt_alloc_align_float(num_elem);
        g->thumb_preview_buf_smoothed_micro = dt_alloc_align_float(num_elem);
        g->thumb_preview_buf_width = width;
        g->thumb_preview_buf_height = height;
        g->luminance_valid = FALSE;
      }

      luminance_pixel = g->thumb_preview_buf_pixel;
      luminance_smoothed = g->thumb_preview_buf_smoothed;
      luminance_smoothed_fine = g->thumb_preview_buf_smoothed_fine;
      luminance_smoothed_micro = g->thumb_preview_buf_smoothed_micro;
      cached = TRUE;
      dt_iop_gui_leave_critical_section(self);
    }
    else
    {
      luminance_pixel = dt_alloc_align_float(num_elem);
      luminance_smoothed = dt_alloc_align_float(num_elem);
      luminance_smoothed_fine = dt_alloc_align_float(num_elem);
      luminance_smoothed_micro = dt_alloc_align_float(num_elem);
    }
  }
  else
  {
    // No interactive editing: allocate local temp buffers
    luminance_pixel = dt_alloc_align_float(num_elem);
    luminance_smoothed = dt_alloc_align_float(num_elem);
    luminance_smoothed_fine = dt_alloc_align_float(num_elem);
    luminance_smoothed_micro = dt_alloc_align_float(num_elem);
  }

  // Check buffer allocation
  if(!luminance_pixel || !luminance_smoothed || !luminance_smoothed_fine || !luminance_smoothed_micro)
  {
    dt_control_log(_("local contrast failed to allocate memory, check your RAM settings"));
    if(!cached)
    {
      dt_free_align(luminance_pixel);
      dt_free_align(luminance_smoothed);
      dt_free_align(luminance_smoothed_fine);
      dt_free_align(luminance_smoothed_micro);
    }
    return;
  }

  // Compute luminance masks
  if(cached)
  {
    if(piece->pipe->type & DT_DEV_PIXELPIPE_FULL)
    {
      dt_hash_t saved_hash;
      hash_set_get(&g->ui_preview_hash, &saved_hash, &self->gui_lock);

      dt_iop_gui_enter_critical_section(self);
      const gboolean luminance_valid = g->luminance_valid;
      dt_iop_gui_leave_critical_section(self);

      if(hash != saved_hash || !luminance_valid)
      {
        compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
        compute_smoothed_luminance_mask(in, luminance_smoothed, width, height, d, d->radius);
        compute_smoothed_luminance_mask(in, luminance_smoothed_fine, width, height, d, d->radius / 2);
        compute_smoothed_luminance_mask(in, luminance_smoothed_micro, width, height, d, d->radius / 4);
        hash_set_get(&hash, &g->ui_preview_hash, &self->gui_lock);
      }
    }
    else if(piece->pipe->type & DT_DEV_PIXELPIPE_PREVIEW)
    {
      dt_hash_t saved_hash;
      hash_set_get(&g->thumb_preview_hash, &saved_hash, &self->gui_lock);

      dt_iop_gui_enter_critical_section(self);
      const gboolean luminance_valid = g->luminance_valid;
      dt_iop_gui_leave_critical_section(self);

      if(saved_hash != hash || !luminance_valid)
      {
        dt_iop_gui_enter_critical_section(self);
        g->thumb_preview_hash = hash;
        compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
        compute_smoothed_luminance_mask(in, luminance_smoothed, width, height, d, d->radius);
        compute_smoothed_luminance_mask(in, luminance_smoothed_fine, width, height, d, d->radius / 2);
        compute_smoothed_luminance_mask(in, luminance_smoothed_micro, width, height, d, d->radius / 4);
        g->luminance_valid = TRUE;
        dt_iop_gui_leave_critical_section(self);
        dt_dev_pixelpipe_cache_invalidate_later(piece->pipe, self->iop_order);
      }
    }
    else
    {
      compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
      compute_smoothed_luminance_mask(in, luminance_smoothed, width, height, d, d->radius);
      compute_smoothed_luminance_mask(in, luminance_smoothed_fine, width, height, d, d->radius / 2);
      compute_smoothed_luminance_mask(in, luminance_smoothed_micro, width, height, d, d->radius / 4);
    }
  }
  else
  {
    compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
    compute_smoothed_luminance_mask(in, luminance_smoothed, width, height, d, d->radius);
    compute_smoothed_luminance_mask(in, luminance_smoothed_fine, width, height, d, d->radius / 2);
    compute_smoothed_luminance_mask(in, luminance_smoothed_micro, width, height, d, d->radius / 4);
  }

  // Display output
  if(self->dev->gui_attached && (piece->pipe->type & DT_DEV_PIXELPIPE_FULL))
  {
    if(g->mask_display)
    {
      display_detail_mask(luminance_pixel, luminance_smoothed, out, width, height);
      piece->pipe->mask_display = DT_DEV_PIXELPIPE_DISPLAY_PASSTHRU;
    }
    else
      apply_local_contrast(in, luminance_pixel, luminance_smoothed,
                           d->fine_scale != 1.0f ? luminance_smoothed_fine : NULL,
                           d->micro_scale != 1.0f ? luminance_smoothed_micro : NULL,
                           out, roi_in, roi_out, d);
  }
  else
  {
    apply_local_contrast(in, luminance_pixel, luminance_smoothed,
                         d->fine_scale != 1.0f ? luminance_smoothed_fine : NULL,
                         d->micro_scale != 1.0f ? luminance_smoothed_micro : NULL,
                         out, roi_in, roi_out, d);
  }

  if(!cached)
  {
    dt_free_align(luminance_pixel);
    dt_free_align(luminance_smoothed);
    dt_free_align(luminance_smoothed_fine);
    dt_free_align(luminance_smoothed_micro);
  }
}


void process(dt_iop_module_t *self,
             dt_dev_pixelpipe_iop_t *piece,
             const void *const restrict ivoid,
             void *const restrict ovoid,
             const dt_iop_roi_t *const roi_in,
             const dt_iop_roi_t *const roi_out)
{
  local_contrast_process(self, piece, ivoid, ovoid, roi_in, roi_out);
}


void modify_roi_in(dt_iop_module_t *self,
                   dt_dev_pixelpipe_iop_t *piece,
                   const dt_iop_roi_t *roi_out,
                   dt_iop_roi_t *roi_in)
{
  dt_iop_local_contrast_rgb_data_t *const d = piece->data;

  // Get the scaled window radius for the box average
  const int max_size = (piece->iwidth > piece->iheight) ? piece->iwidth : piece->iheight;
  const float diameter = d->blending * max_size * roi_in->scale;
  const int radius = (int)((diameter - 1.0f) / 2.0f);
  d->radius = radius;
}


void init_global(dt_iop_module_so_t *self)
{
  dt_iop_local_contrast_rgb_global_data_t *gd = malloc(sizeof(dt_iop_local_contrast_rgb_global_data_t));
  self->data = gd;
}


void cleanup_global(dt_iop_module_so_t *self)
{
  free(self->data);
  self->data = NULL;
}


void commit_params(dt_iop_module_t *self,
                   dt_iop_params_t *p1,
                   dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  const dt_iop_local_contrast_rgb_params_t *p = (dt_iop_local_contrast_rgb_params_t *)p1;
  dt_iop_local_contrast_rgb_data_t *d = piece->data;

  d->method = p->method;
  d->details = p->details;
  d->iterations = p->iterations;
  d->detail_scale = p->detail_scale;
  d->fine_scale = p->fine_scale;
  d->micro_scale = p->micro_scale;
  d->global_scale = p->global_scale;

  // UI blending param is the square root of the actual blending parameter
  // to make it more sensitive to small values that represent the most important value domain.
  // UI parameter is given in percentage of maximum blending value.
  // The actual blending parameter represents the fraction of the largest image dimension.
  d->blending = p->blending * p->blending / 10000.0f;

  // UI guided filter feathering param increases edge taping
  // but actual regularization behaves inversely
  d->feathering = 1.0f / p->feathering;
}


void init_pipe(dt_iop_module_t *self,
               dt_dev_pixelpipe_t *pipe,
               dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = dt_calloc1_align_type(dt_iop_local_contrast_rgb_data_t);
}


void cleanup_pipe(dt_iop_module_t *self,
                  dt_dev_pixelpipe_t *pipe,
                  dt_dev_pixelpipe_iop_t *piece)
{
  dt_free_align(piece->data);
  piece->data = NULL;
}


static void gui_cache_init(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;
  if(g == NULL) return;

  dt_iop_gui_enter_critical_section(self);
  g->ui_preview_hash = DT_INVALID_HASH;
  g->thumb_preview_hash = DT_INVALID_HASH;
  g->mask_display = FALSE;
  g->luminance_valid = FALSE;

  g->full_preview_buf_pixel = NULL;
  g->full_preview_buf_smoothed = NULL;
  g->full_preview_buf_smoothed_fine = NULL;
  g->full_preview_buf_smoothed_micro = NULL;
  g->full_preview_buf_width = 0;
  g->full_preview_buf_height = 0;

  g->thumb_preview_buf_pixel = NULL;
  g->thumb_preview_buf_smoothed = NULL;
  g->thumb_preview_buf_smoothed_fine = NULL;
  g->thumb_preview_buf_smoothed_micro = NULL;
  g->thumb_preview_buf_width = 0;
  g->thumb_preview_buf_height = 0;

  g->pipe_order = 0;
  dt_iop_gui_leave_critical_section(self);
}


static void show_guiding_controls(const dt_iop_module_t *self)
{
  const dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  // All filters need these controls
  gtk_widget_set_visible(g->blending, TRUE);
  gtk_widget_set_visible(g->feathering, TRUE);
  gtk_widget_set_visible(g->iterations, TRUE);
}


void gui_update(dt_iop_module_t *self)
{
  const dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  show_guiding_controls(self);
  invalidate_luminance_cache(self);

  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_luminance_mask), g->mask_display);
}


void gui_changed(dt_iop_module_t *self,
                 GtkWidget *w,
                 void *previous)
{
  const dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  if(w == g->method
     || w == g->blending
     || w == g->feathering
     || w == g->iterations
     || w == g->details)
  {
    invalidate_luminance_cache(self);
  }
}


static void show_luminance_mask_callback(GtkWidget *togglebutton,
                                         GdkEventButton *event,
                                         dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_request_focus(self);

  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->off), TRUE);

  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  // If blend module is displaying mask, don't display here
  if(self->request_mask_display)
  {
    dt_control_log(_("cannot display masks when the blending mask is displayed"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_luminance_mask), FALSE);
    g->mask_display = FALSE;
    return;
  }
  else
    g->mask_display = !g->mask_display;

  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_luminance_mask), g->mask_display);
  dt_iop_refresh_center(self);
}


static void _develop_ui_pipe_started_callback(gpointer instance,
                                              dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;
  if(g == NULL) return;

  if(!self->expanded || !self->enabled)
  {
    dt_iop_gui_enter_critical_section(self);
    g->mask_display = FALSE;
    dt_iop_gui_leave_critical_section(self);
  }

  ++darktable.gui->reset;
  dt_iop_gui_enter_critical_section(self);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_luminance_mask), g->mask_display);
  dt_iop_gui_leave_critical_section(self);
  --darktable.gui->reset;
}


static void _develop_preview_pipe_finished_callback(gpointer instance,
                                                    dt_iop_module_t *self)
{
  const dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;
  if(g == NULL) return;
}


static void _develop_ui_pipe_finished_callback(gpointer instance,
                                               dt_iop_module_t *self)
{
  const dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;
  if(g == NULL) return;
}


void gui_reset(dt_iop_module_t *self)
{
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}


void gui_init(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = IOP_GUI_ALLOC(local_contrast_rgb);

  gui_cache_init(self);

  // Main container
  self->widget = dt_gui_vbox();

  // Detail boost slider
  g->detail_scale = dt_bauhaus_slider_from_params(self, "detail_scale");
  dt_bauhaus_slider_set_soft_range(g->detail_scale, 0.25, 3.0);
  dt_bauhaus_slider_set_digits(g->detail_scale, 2);
  gtk_widget_set_tooltip_text
    (g->detail_scale,
     _("amount of local contrast enhancement\n"
       "1.0 = no change\n"
       "> 1.0 = boost local contrast\n"
       "< 1.0 = reduce local contrast"));

// Fine detail slider
  g->fine_scale = dt_bauhaus_slider_from_params(self, "fine_scale");
  dt_bauhaus_slider_set_soft_range(g->fine_scale, 0.25, 3.0);
  dt_bauhaus_slider_set_digits(g->fine_scale, 2);
  gtk_widget_set_tooltip_text(g->fine_scale, _("amount of fine detail enhancement"));
  dt_gui_box_add(self->widget, g->fine_scale);

  // Micro detail slider
  g->micro_scale = dt_bauhaus_slider_from_params(self, "micro_scale");
  dt_bauhaus_slider_set_soft_range(g->micro_scale, 0.25, 3.0);
  dt_bauhaus_slider_set_digits(g->micro_scale, 2);
  gtk_widget_set_tooltip_text(g->micro_scale, _("amount of micro detail enhancement"));
  dt_gui_box_add(self->widget, g->micro_scale);

  // Global contrast slider
  g->global_scale = dt_bauhaus_slider_from_params(self, "global_scale");
  dt_bauhaus_slider_set_soft_range(g->global_scale, 0.25, 3.0);
  dt_bauhaus_slider_set_digits(g->global_scale, 2);
  gtk_widget_set_tooltip_text
    (g->global_scale,
     _("amount of global contrast enhancement"));

  // Separator
  gtk_widget_set_margin_top(dt_ui_section_label_new(C_("section", "masking")), DT_PIXEL_APPLY_DPI(10));
  dt_gui_box_add(self->widget, dt_ui_section_label_new(C_("section", "masking")));

  g->blending = dt_bauhaus_slider_from_params(self, "blending");
  dt_bauhaus_slider_set_soft_range(g->blending, 0.1, 100.0);
  dt_bauhaus_slider_set_format(g->blending, "%");
  gtk_widget_set_tooltip_text
    (g->blending,
     _("size of the smoothing area as percentage of image size\n"
       "larger = affects broader features\n"
       "smaller = affects finer details"));

  g->feathering = dt_bauhaus_slider_from_params(self, "feathering");
  dt_bauhaus_slider_set_soft_range(g->feathering, 0.1, 50.0);
  gtk_widget_set_tooltip_text
    (g->feathering,
     _("edge sensitivity of the filter\n"
       "higher = better edge preservation\n"
       "lower = smoother transitions, but may lead to halos around edges"));

  // Filter parameters
  g->iterations = dt_bauhaus_slider_from_params(self, "iterations");
  dt_bauhaus_slider_set_soft_max(g->iterations, 5);
  gtk_widget_set_tooltip_text
    (g->iterations,
     _("number of filter passes\n"
       "more iterations = smoother result but slower"));

  // Luminance estimator
  g->method = dt_bauhaus_combobox_from_params(self, "method");
  gtk_widget_set_tooltip_text
    (g->method,
     _("method used to estimate pixel luminance from RGB values\n"
       "choose the one that gives best contrast between details and surroundings"));

  // Detail preservation filter
  g->details = dt_bauhaus_combobox_from_params(self, N_("details"));
  gtk_widget_set_tooltip_text
    (g->details,
     _("edge-aware filter used to smooth the luminance mask\n"
       "'guided filter' is good for general use\n"
       "'EIGF' (exposure-independent guided filter) treats shadows and highlights equally\n"
       "'averaged' variants blend with unfiltered for softer effect"));

  // Display mask toggle
  g->show_luminance_mask = dt_iop_togglebutton_new
    (self, NULL,
     N_("display detail mask"), NULL, G_CALLBACK(show_luminance_mask_callback),
     FALSE, 0, 0, dtgtk_cairo_paint_showmask, NULL);
  dt_gui_add_class(g->show_luminance_mask, "dt_transparent_background");
  dtgtk_togglebutton_set_paint(DTGTK_TOGGLEBUTTON(g->show_luminance_mask),
                               dtgtk_cairo_paint_showmask, 0, NULL);
  dt_gui_add_class(g->show_luminance_mask, "dt_bauhaus_alignment");

  GtkWidget *hbox = dt_gui_hbox(dt_gui_expand(dt_ui_label_new(_("display detail mask"))),
                                g->show_luminance_mask);
  dt_gui_box_add(self->widget, hbox);

  // Connect signals for pipe events
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_PREVIEW_PIPE_FINISHED, _develop_preview_pipe_finished_callback);
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_UI_PIPE_FINISHED, _develop_ui_pipe_finished_callback);
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_HISTORY_CHANGE, _develop_ui_pipe_started_callback);
}


void gui_cleanup(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  dt_free_align(g->thumb_preview_buf_pixel);
  dt_free_align(g->thumb_preview_buf_smoothed);
  dt_free_align(g->thumb_preview_buf_smoothed_fine);
  dt_free_align(g->thumb_preview_buf_smoothed_micro);
  dt_free_align(g->full_preview_buf_pixel);
  dt_free_align(g->full_preview_buf_smoothed);
  dt_free_align(g->full_preview_buf_smoothed_fine);
  dt_free_align(g->full_preview_buf_smoothed_micro);
}


// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
