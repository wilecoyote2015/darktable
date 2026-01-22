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
 * It supports N independent detail scales, each with its own:
 * - Detail boost (local contrast scaling factor)
 * - Feature scale (smoothing diameter)
 * - Edge refinement/feathering
 *
 * The module works by computing luminance masks:
 * 1. A pixel-wise luminance (unblurred) - computed once
 * 2. A smoothed luminance using edge-aware filters - computed per scale
 *
 * The difference between pixel and smoothed luminance represents local detail.
 * Each scale's contribution is calculated independently and summed before
 * applying the final exposure correction.
 *
 * The module should be placed early in the pipe (before color profile)
 * as it operates on scene-linear RGB data.
 *
 ***/

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


DT_MODULE_INTROSPECTION(1, dt_iop_local_contrast_rgb_params_t)


/** Number of independent detail scales */
#define N_SCALES 3

/** Minimum float value to avoid log2(0) */
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
  // Per-scale parameters
  float detail_boost[N_SCALES];   // $MIN: 0.0 $MAX: 500.0 $DEFAULT: 100.0 $DESCRIPTION: "detail boost"
  float feature_scale[N_SCALES];       // $MIN: 0.01 $MAX: 100.0 $DEFAULT: 12.0 $DESCRIPTION: "feature scale"
  float feathering[N_SCALES];     // $MIN: 0.01 $MAX: 10000.0 $DEFAULT: 5.0 $DESCRIPTION: "edges refinement"

  // Shared parameters
  dt_iop_local_contrast_rgb_filter_t details; // $DEFAULT: DT_LC_EIGF $DESCRIPTION: "feature extractor"
  dt_iop_luminance_mask_method_t method;      // $DEFAULT: DT_TONEEQ_NORM_2 $DESCRIPTION: "luminance estimator"
  int iterations;                             // $MIN: 1 $MAX: 20 $DEFAULT: 1 $DESCRIPTION: "filter diffusion"
} dt_iop_local_contrast_rgb_params_t;


/**
 * Per-scale processing data derived from params
 **/
typedef struct dt_iop_local_contrast_rgb_scale_data_t
{
  float detail_boost;
  float feature_scale;
  float feathering;
  int radius;    // derived from feature_scale and image size
} dt_iop_local_contrast_rgb_scale_data_t;


typedef struct dt_iop_local_contrast_rgb_data_t
{
  // Per-scale data
  dt_iop_local_contrast_rgb_scale_data_t scales[N_SCALES];

  // Shared data
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
  // Flags: which scale's mask is displayed (-1 = none)
  int mask_display_scale;

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
  float *thumb_preview_buf_pixel;                    // pixel-wise luminance (no blur)
  float *thumb_preview_buf_smoothed[N_SCALES];       // smoothed luminance per scale
  float *full_preview_buf_pixel;
  float *full_preview_buf_smoothed[N_SCALES];

  // Cache validity
  gboolean luminance_valid;

  // Per-scale GTK widgets
  GtkWidget *detail_boost[N_SCALES];
  GtkWidget *feature_scale[N_SCALES];
  GtkWidget *feathering[N_SCALES];
  GtkWidget *show_mask[N_SCALES];

  // Shared GTK widgets
  GtkWidget *details;
  GtkWidget *iterations;
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
    (self, _("enhance local contrast by boosting fine details while preserving edges\n"
             "supports multiple independent detail scales"),
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
 * Check if any scale is active (detail_boost != 1.0)
 **/
static inline gboolean has_active_scales(const dt_iop_local_contrast_rgb_data_t *const d)
{
  for(int s = 0; s < N_SCALES; s++)
  {
    if(d->scales[s].detail_boost != 1.0f) return TRUE;
  }
  return FALSE;
}


/**
 * Check if a specific scale is active
 **/
static inline gboolean scale_is_active(const float detail_boost)
{
  return detail_boost != 1.0f;
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
 * Compute smoothed luminance mask for a single scale
 **/
__DT_CLONE_TARGETS__
static inline void compute_smoothed_luminance_for_scale(
    const float *const restrict in,
    float *const restrict luminance,
    const size_t width,
    const size_t height,
    const dt_iop_local_contrast_rgb_scale_data_t *const scale_data,
    const dt_iop_local_contrast_rgb_filter_t details,
    const int iterations,
    const dt_iop_luminance_mask_method_t method)
{
  // First compute pixel-wise luminance (no boost)
  luminance_mask(in, luminance, width, height, method, 1.0f, 0.0f, 1.0f);

  // Then apply the smoothing filter
  switch(details)
  {
    case(DT_LC_AVG_GUIDED):
    {
      fast_surface_blur(luminance, width, height,
                        scale_data->radius, scale_data->feathering, iterations,
                        DT_GF_BLENDING_GEOMEAN, 1.0f, 0.0f,
                        exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_GUIDED):
    {
      fast_surface_blur(luminance, width, height,
                        scale_data->radius, scale_data->feathering, iterations,
                        DT_GF_BLENDING_LINEAR, 1.0f, 0.0f,
                        exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_AVG_EIGF):
    {
      fast_eigf_surface_blur(luminance, width, height,
                             scale_data->radius, scale_data->feathering, iterations,
                             DT_GF_BLENDING_GEOMEAN, 1.0f,
                             0.0f, exp2f(-14.0f), 4.0f);
      break;
    }

    case(DT_LC_EIGF):
    {
      fast_eigf_surface_blur(luminance, width, height,
                             scale_data->radius, scale_data->feathering, iterations,
                             DT_GF_BLENDING_LINEAR, 1.0f,
                             0.0f, exp2f(-14.0f), 4.0f);
      break;
    }
  }
}


/**
 * Apply multi-scale local contrast enhancement
 *
 * For each pixel:
 * 1. Sum correction_ev contributions from all active scales
 * 2. Apply single exp2f(total_correction) multiplier
 *
 * This is more efficient than applying each scale sequentially and
 * allows the effects to combine additively in log space.
 **/
__DT_CLONE_TARGETS__
static inline void apply_multiscale_local_contrast(
    const float *const restrict in,
    const float *const restrict luminance_pixel,
    float *const restrict *const restrict luminance_smoothed,
    float *const restrict out,
    const size_t width,
    const size_t height,
    const dt_iop_local_contrast_rgb_data_t *const d)
{
  const size_t npixels = width * height;

  // Unpack scale data for vectorization
  float detail_boosts[N_SCALES] DT_ALIGNED_PIXEL;
  gboolean active[N_SCALES];
  int n_active = 0;

  for(int s = 0; s < N_SCALES; s++)
  {
    detail_boosts[s] = d->scales[s].detail_boost;
    active[s] = scale_is_active(detail_boosts[s]);
    if(active[s]) n_active++;
  }

  // Early exit if no scales are active
  if(n_active == 0)
  {
    dt_iop_image_copy_by_size(out, in, width, height, 4);
    return;
  }

  DT_OMP_FOR()
  for(size_t k = 0; k < npixels; k++)
  {
    const float lum_pixel = fmaxf(luminance_pixel[k], MIN_FLOAT);
    float total_correction_ev = 0.0f;

    // Sum correction contributions from all active scales
    for(int s = 0; s < N_SCALES; s++)
    {
      if(!active[s]) continue;

      const float lum_smoothed = fmaxf(luminance_smoothed[s][k], MIN_FLOAT);

      // Detail in log space (EV): how much brighter/darker is this pixel
      // compared to its local neighborhood at this scale
      const float detail_ev = log2f(lum_pixel / lum_smoothed);

      // Scale the detail: detail_boost = 1.0 means no change
      // > 1.0 boosts local contrast, < 1.0 reduces it
      const float scaled_detail_ev = detail_boosts[s] * detail_ev;

      // The correction is the difference between scaled and original detail
      total_correction_ev += scaled_detail_ev - detail_ev;
    }

    // Apply combined correction in linear space
    const float multiplier = exp2f(total_correction_ev);

    for_each_channel(c)
      out[4 * k + c] = in[4 * k + c] * multiplier;
  }
}


/**
 * Display the detail mask for a specific scale
 * Output is a grayscale image normalized to [0, 1] where:
 * - 0.5 = no local detail (pixel matches neighborhood)
 * - < 0.5 = pixel darker than neighborhood
 * - > 0.5 = pixel brighter than neighborhood
 **/
__DT_CLONE_TARGETS__
static inline void display_detail_mask_for_scale(
    const float *const restrict luminance_pixel,
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
    out[4 * k + 3] = 1.0f;
  }
}


/**
 * Allocate or reallocate smoothed buffers for all scales
 **/
static inline void alloc_smoothed_buffers(float **buffers,
                                          const size_t num_elem)
{
  for(int s = 0; s < N_SCALES; s++)
  {
    dt_free_align(buffers[s]);
    buffers[s] = dt_alloc_align_float(num_elem);
  }
}


/**
 * Free smoothed buffers for all scales
 **/
static inline void free_smoothed_buffers(float **buffers)
{
  for(int s = 0; s < N_SCALES; s++)
  {
    dt_free_align(buffers[s]);
    buffers[s] = NULL;
  }
}


/**
 * Check if all smoothed buffers are allocated
 **/
static inline gboolean smoothed_buffers_valid(float **buffers)
{
  for(int s = 0; s < N_SCALES; s++)
  {
    if(!buffers[s]) return FALSE;
  }
  return TRUE;
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
  float *luminance_smoothed[N_SCALES] = { NULL };

  const size_t width = roi_in->width;
  const size_t height = roi_in->height;
  const size_t num_elem = width * height;

  // Get the hash of the upstream pipe to track changes
  const dt_hash_t hash = dt_dev_pixelpipe_piece_hash(piece, roi_out, TRUE);

  // Sanity checks
  if(width < 1 || height < 1) return;
  if(roi_in->width < roi_out->width || roi_in->height < roi_out->height) return;
  if(piece->colors != 4) return;

  // Fast path: if no scales are active, just copy
  if(!has_active_scales(d))
  {
    dt_iop_image_copy_by_size(out, in, width, height, 4);
    return;
  }

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
        g->full_preview_buf_pixel = dt_alloc_align_float(num_elem);
        alloc_smoothed_buffers(g->full_preview_buf_smoothed, num_elem);
        g->full_preview_buf_width = width;
        g->full_preview_buf_height = height;
      }

      luminance_pixel = g->full_preview_buf_pixel;
      for(int s = 0; s < N_SCALES; s++)
        luminance_smoothed[s] = g->full_preview_buf_smoothed[s];
      cached = TRUE;
    }
    else if(piece->pipe->type & DT_DEV_PIXELPIPE_PREVIEW)
    {
      dt_iop_gui_enter_critical_section(self);
      if(g->thumb_preview_buf_width != width || g->thumb_preview_buf_height != height)
      {
        dt_free_align(g->thumb_preview_buf_pixel);
        g->thumb_preview_buf_pixel = dt_alloc_align_float(num_elem);
        alloc_smoothed_buffers(g->thumb_preview_buf_smoothed, num_elem);
        g->thumb_preview_buf_width = width;
        g->thumb_preview_buf_height = height;
        g->luminance_valid = FALSE;
      }

      luminance_pixel = g->thumb_preview_buf_pixel;
      for(int s = 0; s < N_SCALES; s++)
        luminance_smoothed[s] = g->thumb_preview_buf_smoothed[s];
      cached = TRUE;
      dt_iop_gui_leave_critical_section(self);
    }
    else
    {
      luminance_pixel = dt_alloc_align_float(num_elem);
      alloc_smoothed_buffers(luminance_smoothed, num_elem);
    }
  }
  else
  {
    // No interactive editing: allocate local temp buffers
    luminance_pixel = dt_alloc_align_float(num_elem);
    alloc_smoothed_buffers(luminance_smoothed, num_elem);
  }

  // Check buffer allocation
  if(!luminance_pixel || !smoothed_buffers_valid(luminance_smoothed))
  {
    dt_control_log(_("local contrast failed to allocate memory, check your RAM settings"));
    if(!cached)
    {
      dt_free_align(luminance_pixel);
      free_smoothed_buffers(luminance_smoothed);
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
        // Compute pixel luminance once
        compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);

        // Compute smoothed luminance for each active scale
        for(int s = 0; s < N_SCALES; s++)
        {
          if(scale_is_active(d->scales[s].detail_boost))
          {
            compute_smoothed_luminance_for_scale(in, luminance_smoothed[s],
                                                 width, height,
                                                 &d->scales[s],
                                                 d->details, d->iterations,
                                                 d->method);
          }
        }
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

        // Compute pixel luminance once
        compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);

        // Compute smoothed luminance for each active scale
        for(int s = 0; s < N_SCALES; s++)
        {
          if(scale_is_active(d->scales[s].detail_boost))
          {
            compute_smoothed_luminance_for_scale(in, luminance_smoothed[s],
                                                 width, height,
                                                 &d->scales[s],
                                                 d->details, d->iterations,
                                                 d->method);
          }
        }

        g->luminance_valid = TRUE;
        dt_iop_gui_leave_critical_section(self);
        dt_dev_pixelpipe_cache_invalidate_later(piece->pipe, self->iop_order);
      }
    }
    else
    {
      // Non-cached pipe: compute everything
      compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
      for(int s = 0; s < N_SCALES; s++)
      {
        if(scale_is_active(d->scales[s].detail_boost))
        {
          compute_smoothed_luminance_for_scale(in, luminance_smoothed[s],
                                               width, height,
                                               &d->scales[s],
                                               d->details, d->iterations,
                                               d->method);
        }
      }
    }
  }
  else
  {
    // Non-GUI: compute everything
    compute_pixel_luminance_mask(in, luminance_pixel, width, height, d->method);
    for(int s = 0; s < N_SCALES; s++)
    {
      if(scale_is_active(d->scales[s].detail_boost))
      {
        compute_smoothed_luminance_for_scale(in, luminance_smoothed[s],
                                             width, height,
                                             &d->scales[s],
                                             d->details, d->iterations,
                                             d->method);
      }
    }
  }

  // Display output
  if(self->dev->gui_attached && (piece->pipe->type & DT_DEV_PIXELPIPE_FULL))
  {
    const int display_scale = g->mask_display_scale;
    if(display_scale >= 0 && display_scale < N_SCALES
       && scale_is_active(d->scales[display_scale].detail_boost)
       && luminance_smoothed[display_scale])
    {
      display_detail_mask_for_scale(luminance_pixel, luminance_smoothed[display_scale],
                                    out, width, height);
      piece->pipe->mask_display = DT_DEV_PIXELPIPE_DISPLAY_PASSTHRU;
    }
    else
    {
      apply_multiscale_local_contrast(in, luminance_pixel, luminance_smoothed,
                                      out, width, height, d);
    }
  }
  else
  {
    apply_multiscale_local_contrast(in, luminance_pixel, luminance_smoothed,
                                    out, width, height, d);
  }

  if(!cached)
  {
    dt_free_align(luminance_pixel);
    free_smoothed_buffers(luminance_smoothed);
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

  // Get the scaled window radius for each scale
  const int max_size = (piece->iwidth > piece->iheight) ? piece->iwidth : piece->iheight;

  for(int s = 0; s < N_SCALES; s++)
  {
    const float diameter = d->scales[s].feature_scale * max_size * roi_in->scale;
    const int radius = (int)((diameter - 1.0f) / 2.0f);
    d->scales[s].radius = radius;
  }
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


void reload_defaults(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_params_t *d = self->default_params;

  // Set 12% scale to have a visible effect by default
  // because this resembles a clarity-like effect
  // Other scales remain at 100% (no effect)
  for(int s = 0; s < N_SCALES; s++)
    if(s != 1)
      d->detail_boost[s] = 100.0f;
    else
      d->detail_boost[s] = 150.0f;

  d->feature_scale[0] = 4.0f;
  d->feature_scale[1] = 12.0f;
  d->feature_scale[2] = 25.0f;
}


void commit_params(dt_iop_module_t *self,
                   dt_iop_params_t *p1,
                   dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  const dt_iop_local_contrast_rgb_params_t *p = (dt_iop_local_contrast_rgb_params_t *)p1;
  dt_iop_local_contrast_rgb_data_t *d = piece->data;

  // Copy shared params
  d->method = p->method;
  d->details = p->details;
  d->iterations = p->iterations;

  // Copy per-scale params and compute derived values
  for(int s = 0; s < N_SCALES; s++)
  {
    // UI parameter is given in percentage of detail strength, where 100% means no change 
    // and 0% means that detail is removed (multiplier 0). Internal math is a multiplier of relative detail EV
    // so that 100% means no change, 200% means double the detail, 50% means half the detail, 
    d->scales[s].detail_boost = p->detail_boost[s] / 100.0f;

    // UI feature_scale param is the square root of the actual feature_scale parameter
    // to make it more sensitive to small values that represent the most important value domain.
    // UI parameter is given in percentage of maximum feature_scale value.
    // The actual feature_scale parameter represents the fraction of the largest image dimension.
    d->scales[s].feature_scale = p->feature_scale[s] * p->feature_scale[s] / 10000.0f;

    // UI guided filter feathering param increases edge taping
    // but actual regularization behaves inversely
    d->scales[s].feathering = 1.0f / p->feathering[s];
  }
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
  g->mask_display_scale = -1;  // no mask displayed
  g->luminance_valid = FALSE;

  g->full_preview_buf_pixel = NULL;
  for(int s = 0; s < N_SCALES; s++)
    g->full_preview_buf_smoothed[s] = NULL;
  g->full_preview_buf_width = 0;
  g->full_preview_buf_height = 0;

  g->thumb_preview_buf_pixel = NULL;
  for(int s = 0; s < N_SCALES; s++)
    g->thumb_preview_buf_smoothed[s] = NULL;
  g->thumb_preview_buf_width = 0;
  g->thumb_preview_buf_height = 0;

  g->pipe_order = 0;
  dt_iop_gui_leave_critical_section(self);
}


void gui_update(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  invalidate_luminance_cache(self);

  // Update mask toggle buttons
  for(int s = 0; s < N_SCALES; s++)
  {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_mask[s]),
                                 g->mask_display_scale == s);
  }
}


void gui_changed(dt_iop_module_t *self,
                 GtkWidget *w,
                 void *previous)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  // Check if any masking-related widget changed
  gboolean invalidate = FALSE;

  // Check shared widgets
  if(w == g->details || w == g->iterations)
  {
    invalidate = TRUE;
  }

  // Check per-scale widgets
  for(int s = 0; s < N_SCALES && !invalidate; s++)
  {
    if(w == g->feature_scale[s] || w == g->feathering[s])
    {
      invalidate = TRUE;
    }
  }

  if(invalidate)
  {
    invalidate_luminance_cache(self);
  }
}


/**
 * Callback for mask display toggle buttons.
 * Ensures only one mask can be displayed at a time.
 **/
static void show_mask_callback(GtkWidget *togglebutton,
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
    dt_control_log(_("cannot display masks when the feature_scale mask is displayed"));
    for(int s = 0; s < N_SCALES; s++)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_mask[s]), FALSE);
    g->mask_display_scale = -1;
    return;
  }

  // Find which scale's button was clicked
  int clicked_scale = -1;
  for(int s = 0; s < N_SCALES; s++)
  {
    if(togglebutton == g->show_mask[s])
    {
      clicked_scale = s;
      break;
    }
  }

  // Toggle: if same scale was active, turn off; otherwise switch to new scale
  if(g->mask_display_scale == clicked_scale)
  {
    g->mask_display_scale = -1;  // turn off
  }
  else
  {
    g->mask_display_scale = clicked_scale;  // switch to new
  }

  // Update all toggle buttons
  for(int s = 0; s < N_SCALES; s++)
  {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_mask[s]),
                                 g->mask_display_scale == s);
  }

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
    g->mask_display_scale = -1;
    dt_iop_gui_leave_critical_section(self);
  }

  ++darktable.gui->reset;
  dt_iop_gui_enter_critical_section(self);
  for(int s = 0; s < N_SCALES; s++)
  {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->show_mask[s]),
                                 g->mask_display_scale == s);
  }
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


/**
 * Create GUI widgets for a single scale section
 **/
static void create_scale_section(dt_iop_module_t *self,
                                 dt_iop_local_contrast_rgb_gui_data_t *g,
                                 const int scale_idx)
{
  // Section label
  char section_name[64];
  snprintf(section_name, sizeof(section_name), _("detail scale %d"), scale_idx + 1);

  gtk_widget_set_margin_top(dt_ui_section_label_new(section_name), DT_PIXEL_APPLY_DPI(10));
  dt_gui_box_add(self->widget, dt_ui_section_label_new(section_name));

  // Detail boost slider
  char param_name[64];
  snprintf(param_name, sizeof(param_name), "detail_boost[%d]", scale_idx);
  g->detail_boost[scale_idx] = dt_bauhaus_slider_from_params(self, param_name);
  dt_bauhaus_slider_set_soft_range(g->detail_boost[scale_idx], 0.0, 300.0);
  dt_bauhaus_slider_set_format(g->detail_boost[scale_idx], "%");
  dt_bauhaus_slider_set_digits(g->detail_boost[scale_idx], 2);
  dt_bauhaus_widget_set_label(g->detail_boost[scale_idx], NULL, _("detail strength"));
  gtk_widget_set_tooltip_text
    (g->detail_boost[scale_idx],
     _("amount of local contrast for this scale\n"
       "100% = no change\n"
       "> 100% = increase local contrast\n"
       "< 100% = decrease local contrast"));

  // Feature scale slider
  snprintf(param_name, sizeof(param_name), "feature_scale[%d]", scale_idx);
  g->feature_scale[scale_idx] = dt_bauhaus_slider_from_params(self, param_name);
  dt_bauhaus_slider_set_soft_range(g->feature_scale[scale_idx], 0.1, 100.0);
  dt_bauhaus_slider_set_format(g->feature_scale[scale_idx], "%");
  dt_bauhaus_widget_set_label(g->feature_scale[scale_idx], NULL, _("feature scale"));
  gtk_widget_set_tooltip_text
    (g->feature_scale[scale_idx],
     _("size of the smoothing area as percentage of image size\n"
       "larger = affects broader features\n"
       "smaller = affects finer details"));

  // Edge refinement slider
  snprintf(param_name, sizeof(param_name), "feathering[%d]", scale_idx);
  g->feathering[scale_idx] = dt_bauhaus_slider_from_params(self, param_name);
  dt_bauhaus_slider_set_soft_range(g->feathering[scale_idx], 0.1, 50.0);
  dt_bauhaus_widget_set_label(g->feathering[scale_idx], NULL, _("edges refinement"));
  gtk_widget_set_tooltip_text
    (g->feathering[scale_idx],
     _("edge sensitivity of the filter\n"
       "higher = better edge preservation\n"
       "lower = smoother transitions, but may lead to halos around edges"));

  // Mask display toggle
  g->show_mask[scale_idx] = dt_iop_togglebutton_new
    (self, NULL,
     N_("display detail mask"), NULL, G_CALLBACK(show_mask_callback),
     FALSE, 0, 0, dtgtk_cairo_paint_showmask, NULL);
  dt_gui_add_class(g->show_mask[scale_idx], "dt_transparent_background");
  dtgtk_togglebutton_set_paint(DTGTK_TOGGLEBUTTON(g->show_mask[scale_idx]),
                               dtgtk_cairo_paint_showmask, 0, NULL);
  dt_gui_add_class(g->show_mask[scale_idx], "dt_bauhaus_alignment");

  GtkWidget *hbox = dt_gui_hbox(dt_gui_expand(dt_ui_label_new(_("display detail mask"))),
                                g->show_mask[scale_idx]);
  dt_gui_box_add(self->widget, hbox);
}


void gui_init(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = IOP_GUI_ALLOC(local_contrast_rgb);

  gui_cache_init(self);

  // Main container
  self->widget = dt_gui_vbox();

  // Create GUI for each scale
  for(int s = 0; s < N_SCALES; s++)
  {
    create_scale_section(self, g, s);
  }

  // Masking section (shared parameters)
  gtk_widget_set_margin_top(dt_ui_section_label_new(C_("section", "masking")), DT_PIXEL_APPLY_DPI(10));
  dt_gui_box_add(self->widget, dt_ui_section_label_new(C_("section", "masking")));

  // Feature extractor
  g->details = dt_bauhaus_combobox_from_params(self, N_("details"));
  gtk_widget_set_tooltip_text
    (g->details,
     _("edge-aware filter used to smooth the luminance mask\n"
       "'guided filter' is good for general use\n"
       "'EIGF' (exposure-independent guided filter) treats shadows and highlights equally\n"
       "'averaged' variants blend with unfiltered for softer effect"));

  // Filter diffusion
  g->iterations = dt_bauhaus_slider_from_params(self, "iterations");
  dt_bauhaus_slider_set_soft_max(g->iterations, 5);
  gtk_widget_set_tooltip_text
    (g->iterations,
     _("number of filter passes\n"
       "more iterations = smoother result but slower"));

  // Connect signals for pipe events
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_PREVIEW_PIPE_FINISHED, _develop_preview_pipe_finished_callback);
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_UI_PIPE_FINISHED, _develop_ui_pipe_finished_callback);
  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_DEVELOP_HISTORY_CHANGE, _develop_ui_pipe_started_callback);
}


void gui_cleanup(dt_iop_module_t *self)
{
  dt_iop_local_contrast_rgb_gui_data_t *g = self->gui_data;

  dt_free_align(g->thumb_preview_buf_pixel);
  free_smoothed_buffers(g->thumb_preview_buf_smoothed);
  dt_free_align(g->full_preview_buf_pixel);
  free_smoothed_buffers(g->full_preview_buf_smoothed);
}


// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
