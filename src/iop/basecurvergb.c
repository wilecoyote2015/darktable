/*
    This file is part of darktable,
    Copyright (C) 2010-2025 darktable developers.

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

#include "bauhaus/bauhaus.h"
#include "common/colorspaces_inline_conversions.h"
#include "common/chromatic_adaptation.h"
#include "common/debug.h"
#include "common/imagebuf.h"
#include "common/math.h"
#include "common/opencl.h"
#include "common/rgb_norms.h"
#include "control/control.h"
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
#include "develop/imageop_gui.h"
#include "develop/tiling.h"
#include "dtgtk/drawingarea.h"
#include "gui/draw.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include "gui/accelerators.h"
#include "iop/iop_api.h"

#include <regex.h>
#include <assert.h>
#include <gtk/gtk.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#define DT_GUI_CURVE_EDITOR_INSET DT_PIXEL_APPLY_DPI(5)
#define DT_IOP_TONECURVE_RES 256
#define MAXNODES 20


DT_MODULE_INTROSPECTION(1, dt_iop_basecurvergb_params_t)

typedef struct dt_iop_basecurvergb_node_t
{
  float x; // $MIN: 0.0 $MAX: 1.0
  float y; // $MIN: 0.0 $MAX: 1.0
} dt_iop_basecurvergb_node_t;

typedef struct dt_iop_basecurvergb_params_t
{
  // three curves (c, ., .) with max number of nodes
  // the other two are reserved, maybe we'll have cam rgb at some point.
  dt_iop_basecurvergb_node_t basecurve[3][MAXNODES];
  int basecurvergb_nodes[3]; // $MIN: 0 $MAX: MAXNODES $DEFAULT: 0
  int basecurvergb_type[3];  // $MIN: 0 $MAX: MONOTONE_HERMITE $DEFAULT: MONOTONE_HERMITE
  float preserve_hue;    // $MIN: 0.0 $MAX: 1.0 $DEFAULT: 0.0 $DESCRIPTION: "preserve hue after application of base curve"
  float preserve_highlight_saturation; // $MIN: 0.0 $MAX: 1.0 $DEFAULT: 0.0 $DESCRIPTION: "preserve highlight saturation after application of base curve"
  float source_white;    // $MIN: -4 $MAX: 4.0 $DEFAULT: 0.0 $DESCRIPTION: "exposure shift before curve is applied"
} dt_iop_basecurvergb_params_t;

int legacy_params(dt_iop_module_t *self,
                  const void *const old_params,
                  const int old_version,
                  void **new_params,
                  int32_t *new_params_size,
                  int *new_version)
{

  return 1;
}

typedef struct dt_iop_basecurvergb_gui_data_t
{
  dt_draw_curve_t *minmax_curve; // curve for gui to draw
  int minmax_curve_type, minmax_curve_nodes;
  GtkBox *hbox;
  GtkDrawingArea *area;
  GtkWidget *preserve_hue;
  GtkWidget *preserve_highlight_saturation;
  GtkWidget *source_white;
  double mouse_x, mouse_y;
  int selected;
  double selected_offset, selected_y, selected_min, selected_max;
  float draw_xs[DT_IOP_TONECURVE_RES], draw_ys[DT_IOP_TONECURVE_RES];
  float draw_min_xs[DT_IOP_TONECURVE_RES], draw_min_ys[DT_IOP_TONECURVE_RES];
  float draw_max_xs[DT_IOP_TONECURVE_RES], draw_max_ys[DT_IOP_TONECURVE_RES];
  float loglogscale;
  GtkWidget *logbase;
} dt_iop_basecurvergb_gui_data_t;

typedef struct basecurvergb_preset_t
{
  const char *name;
  const char *maker;
  const char *model;
  int iso_min;
  float iso_max;
  dt_iop_basecurvergb_params_t params;
  int filter;
} basecurvergb_preset_t;

#define m MONOTONE_HERMITE

static const basecurvergb_preset_t basecurvergb_camera_presets[] = {
  // copy paste your measured basecurve line at the top here, like so (note the exif data and the last 1):
  // clang-format off

  // nikon d750 by Edouard Gomez
  {"Nikon D750", "NIKON CORPORATION", "NIKON D750", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.018124, 0.026126}, {0.143357, 0.370145}, {0.330116, 0.730507}, {0.457952, 0.853462}, {0.734950, 0.965061}, {0.904758, 0.985699}, {1.000000, 1.000000}}}, {8}, {m}, 1, 0, 0}, 1},
  // contributed by Stefan Kauerauf
  {"Nikon D5100", "NIKON CORPORATION", "NIKON D5100", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.001113, 0.000506}, {0.002842, 0.001338}, {0.005461, 0.002470}, {0.011381, 0.006099}, {0.013303, 0.007758}, {0.034638, 0.041119}, {0.044441, 0.063882}, {0.070338, 0.139639}, {0.096068, 0.210915}, {0.137693, 0.310295}, {0.206041, 0.432674}, {0.255508, 0.504447}, {0.302770, 0.569576}, {0.425625, 0.726755}, {0.554526, 0.839541}, {0.621216, 0.882839}, {0.702662, 0.927072}, {0.897426, 0.990984}, {1.000000, 1.000000}}}, {20}, {m}, 1, 0, 0}, 1},
  // nikon d7000 by Edouard Gomez
  {"Nikon D7000", "NIKON CORPORATION", "NIKON D7000", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.001943, 0.003040}, {0.019814, 0.028810}, {0.080784, 0.210476}, {0.145700, 0.383873}, {0.295961, 0.654041}, {0.651915, 0.952819}, {1.000000, 1.000000}}}, {8}, {m}, 1, 0, 0}, 1},
  // nikon d7200 standard by Ralf Brown (firmware 1.00)
  {"Nikon D7200", "NIKON CORPORATION", "NIKON D7200", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.001604, 0.001334}, {0.007401, 0.005237}, {0.009474, 0.006890}, {0.017348, 0.017176}, {0.032782, 0.044336}, {0.048033, 0.086548}, {0.075803, 0.168331}, {0.109539, 0.273539}, {0.137373, 0.364645}, {0.231651, 0.597511}, {0.323797, 0.736475}, {0.383796, 0.805797}, {0.462284, 0.872247}, {0.549844, 0.918328}, {0.678855, 0.962361}, {0.817445, 0.990406}, {1.000000, 1.000000}}}, {18}, {m}, 1, 0, 0}, 1},
  // nikon d7500 by Anders Bennehag (firmware C 1.00, LD 2.016)
  {"NIKON D7500", "NIKON CORPORATION", "NIKON D7500", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.000892, 0.001062}, {0.002280, 0.001768}, {0.013983, 0.011368}, {0.032597, 0.044700}, {0.050065, 0.097131}, {0.084129, 0.219954}, {0.120975, 0.336806}, {0.170730, 0.473752}, {0.258677, 0.647113}, {0.409997, 0.827417}, {0.499979, 0.889468}, {0.615564, 0.941960}, {0.665272, 0.957736}, {0.832126, 0.991968}, {1.000000, 1.000000}}}, {16}, {m}, 1, 0, 0}, 1},
  // sony rx100m2 by Günther R.
  { "Sony DSC-RX100M2", "SONY", "DSC-RX100M2", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.015106, 0.008116 }, { 0.070077, 0.093725 }, { 0.107484, 0.170723 }, { 0.191528, 0.341093 }, { 0.257996, 0.458453 }, { 0.305381, 0.537267 }, { 0.326367, 0.569257 }, { 0.448067, 0.723742 }, { 0.509627, 0.777966 }, { 0.676751, 0.898797 }, { 1.000000, 1.000000 } } }, { 12 }, { m }, 1}, 1 },
  // contributed by matthias bodenbinder
  { "Canon EOS 6D", "Canon", "Canon EOS 6D", 0, FLT_MAX, { { { { 0.000000, 0.002917 }, { 0.000751, 0.001716 }, { 0.006011, 0.004438 }, { 0.020286, 0.021725 }, { 0.048084, 0.085918 }, { 0.093914, 0.233804 }, { 0.162284, 0.431375 }, { 0.257701, 0.629218 }, { 0.384673, 0.800332 }, { 0.547709, 0.917761 }, { 0.751315, 0.988132 }, { 1.000000, 0.999943 } } }, { 12 }, { m }, 1}, 1 },
  // contributed by Dan Torop
  { "Fujifilm X100S", "Fujifilm", "X100S", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.009145, 0.007905 }, { 0.026570, 0.032201 }, { 0.131526, 0.289717 }, { 0.175858, 0.395263 }, { 0.350981, 0.696899 }, { 0.614997, 0.959451 }, { 1.000000, 1.000000 } } }, { 8 }, { m }, 1}, 1 },
  { "Fujifilm X100T", "Fujifilm", "X100T", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.009145, 0.007905 }, { 0.026570, 0.032201 }, { 0.131526, 0.289717 }, { 0.175858, 0.395263 }, { 0.350981, 0.696899 }, { 0.614997, 0.959451 }, { 1.000000, 1.000000 } } }, { 8 }, { m }, 1}, 1 },
  // contributed by Johannes Hanika
  { "Canon EOS 5D Mark II", "Canon", "Canon EOS 5D Mark II", 0, FLT_MAX, { { { { 0.000000, 0.000366 }, { 0.006560, 0.003504 }, { 0.027310, 0.029834 }, { 0.045915, 0.070230 }, { 0.206554, 0.539895 }, { 0.442337, 0.872409 }, { 0.673263, 0.971703 }, { 1.000000, 0.999832 } } }, { 8 }, { m }, 1}, 1 },
  // contributed by chrik5
  { "Pentax K-5", "Pentax", "Pentax K-5", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.004754, 0.002208 }, { 0.009529, 0.004214 }, { 0.023713, 0.013508 }, { 0.031866, 0.020352 }, { 0.046734, 0.034063 }, { 0.059989, 0.052413 }, { 0.088415, 0.096030 }, { 0.136610, 0.190629 }, { 0.174480, 0.256484 }, { 0.205192, 0.307430 }, { 0.228896, 0.348447 }, { 0.286411, 0.428680 }, { 0.355314, 0.513527 }, { 0.440014, 0.607651 }, { 0.567096, 0.732791 }, { 0.620597, 0.775968 }, { 0.760355, 0.881828 }, { 0.875139, 0.960682 }, { 1.000000, 1.000000 } } }, { 20 }, { m }, 1}, 1 },
  // contributed by Togan Muftuoglu - ed: slope is too aggressive on shadows
  //{ "Nikon D90", "NIKON", "D90", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.015520, 0.012248 }, { 0.097950, 0.251013 }, { 0.301515, 0.621951 }, { 0.415513, 0.771384 }, { 0.547326, 0.843079 }, { 0.819769, 0.956678 }, { 1.000000, 1.000000 } } }, { 8 }, { m }, 1, 0, 0 }, 0, 1 },
  // contributed by Edouard Gomez
  {"Nikon D90", "NIKON CORPORATION", "NIKON D90", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.011702, 0.012659}, {0.122918, 0.289973}, {0.153642, 0.342731}, {0.246855, 0.510114}, {0.448958, 0.733820}, {0.666759, 0.894290}, {1.000000, 1.000000}}}, {8}, {m}, 1, 0, 0}, 1},
  // contributed by Pascal Obry
  { "Nikon D800", "NIKON", "NIKON D800", 0, FLT_MAX, { { { { 0.000000, 0.000000 }, { 0.001773, 0.001936 }, { 0.009671, 0.009693 }, { 0.016754, 0.020617 }, { 0.024884, 0.037309 }, { 0.048174, 0.107768 }, { 0.056932, 0.139532 }, { 0.085504, 0.233303 }, { 0.130378, 0.349747 }, { 0.155476, 0.405445 }, { 0.175245, 0.445918 }, { 0.217657, 0.516873 }, { 0.308475, 0.668608 }, { 0.375381, 0.754058 }, { 0.459858, 0.839909 }, { 0.509567, 0.881543 }, { 0.654394, 0.960877 }, { 0.783380, 0.999161 }, { 0.859310, 1.000000 }, { 1.000000, 1.000000 } } }, { 20 }, { m }, 1}, 1 },
  // contributed by Lukas Schrangl
  {"Olympus OM-D E-M10 II", "OLYMPUS CORPORATION    ", "E-M10MarkII     ", 0, FLT_MAX, {{{{0.000000, 0.000000}, {0.005707, 0.004764}, {0.018944, 0.024456}, {0.054501, 0.129992}, {0.075665, 0.211873}, {0.119641, 0.365771}, {0.173148, 0.532024}, {0.247979, 0.668989}, {0.357597, 0.780138}, {0.459003, 0.839829}, {0.626844, 0.904426}, {0.769425, 0.948541}, {0.820429, 0.964715}, {1.000000, 1.000000}}}, {14}, {m}, 1, 0, 0}, 1},
  // clang-format on
};
static const int basecurvergb_camera_presets_cnt = sizeof(basecurvergb_camera_presets) / sizeof(basecurvergb_preset_t);

static const basecurvergb_preset_t basecurvergb_presets[] = {
  // clang-format off
  // smoother cubic spline curve
  { N_("cubic spline"),             "", "",                      0, FLT_MAX, { { { { 0.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } , { 0.0, 0.0}, { 0.0, 0.0} } }            , { 2 }, { CUBIC_SPLINE } }, 0 },
  { N_("neutral"),                  "", "",                      0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.005000, 0.002500 }, { 0.150000, 0.300000 }, { 0.400000, 0.700000 }, { 0.750000, 0.950000 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 1 },
  { N_("canon eos like"),           "Canon", "",                 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.028226, 0.029677 }, { 0.120968, 0.232258 }, { 0.459677, 0.747581 }, { 0.858871, 0.967742 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("canon eos like alternate"), "Canon", "EOS 5D Mark%",     0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.026210, 0.029677 }, { 0.108871, 0.232258 }, { 0.350806, 0.747581 }, { 0.669355, 0.967742 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("nikon like"),               "NIKON", "",                 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.036290, 0.036532 }, { 0.120968, 0.228226 }, { 0.459677, 0.759678 }, { 0.858871, 0.983468 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("nikon like alternate"),     "NIKON", "%D____%",          0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.012097, 0.007322 }, { 0.072581, 0.130742 }, { 0.310484, 0.729291 }, { 0.611321, 0.951613 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("sony alpha like"),          "SONY", "",                  0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.031949, 0.036532 }, { 0.105431, 0.228226 }, { 0.434505, 0.759678 }, { 0.855738, 0.983468 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("pentax like"),              "PENTAX", "",                0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.032258, 0.024596 }, { 0.120968, 0.166419 }, { 0.205645, 0.328527 }, { 0.604839, 0.790171 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("ricoh like"),               "RICOH", "",                 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.032259, 0.024596 }, { 0.120968, 0.166419 }, { 0.205645, 0.328527 }, { 0.604839, 0.790171 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("olympus like"),             "OLYMPUS", "",               0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.033962, 0.028226 }, { 0.249057, 0.439516 }, { 0.501887, 0.798387 }, { 0.750943, 0.955645 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("olympus like alternate"),   "OLYMPUS", "E-M%",           0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.012097, 0.010322 }, { 0.072581, 0.167742 }, { 0.310484, 0.711291 }, { 0.645161, 0.956855 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("panasonic like"),           "Panasonic", "",             0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.036290, 0.024596 }, { 0.120968, 0.166419 }, { 0.205645, 0.328527 }, { 0.604839, 0.790171 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("leica like"),               "Leica", "",                 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.036291, 0.024596 }, { 0.120968, 0.166419 }, { 0.205645, 0.328527 }, { 0.604839, 0.790171 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("kodak easyshare like"),     "EASTMAN KODAK COMPANY", "", 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.044355, 0.020967 }, { 0.133065, 0.154322 }, { 0.209677, 0.300301 }, { 0.572581, 0.753477 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("konica minolta like"),      "MINOLTA", "",               0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.020161, 0.010322 }, { 0.112903, 0.167742 }, { 0.500000, 0.711291 }, { 0.899194, 0.956855 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("samsung like"),             "SAMSUNG", "",               0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.040323, 0.029677 }, { 0.133065, 0.232258 }, { 0.447581, 0.747581 }, { 0.842742, 0.967742 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("fujifilm like"),            "FUJIFILM", "",              0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.028226, 0.029677 }, { 0.104839, 0.232258 }, { 0.387097, 0.747581 }, { 0.754032, 0.967742 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  { N_("nokia like"),               "Nokia", "",                 0, FLT_MAX, { { { { 0.0, 0.0 }, { 0.041825, 0.020161 }, { 0.117871, 0.153226 }, { 0.319392, 0.500000 }, { 0.638783, 0.842742 }, { 1.000000, 1.000000 } } }, { 6 }, { m }, 1, 0, 0 }, 0 },
  // clang-format on
};
#undef m
static const int basecurvergb_presets_cnt = sizeof(basecurvergb_presets) / sizeof(basecurvergb_preset_t);

typedef struct dt_iop_basecurvergb_data_t
{
  dt_draw_curve_t *curve; // curve for pixelpipe piece and pixel processing
  int basecurvergb_type;
  int basecurvergb_nodes;
  float table[0x10000];      // precomputed look-up table for tone curve
  float unbounded_coeffs[3]; // approximation for extrapolation
  float preserve_hue;
  float preserve_highlight_saturation;
  float source_white;
} dt_iop_basecurvergb_data_t;

typedef struct dt_iop_basecurvergb_global_data_t
{
  int kernel_basecurvergb_lut;
} dt_iop_basecurvergb_global_data_t;



const char *name()
{
  return _("base curve rgb");
}

const char **description(dt_iop_module_t *self)
{
  return dt_iop_set_description
    (self,
     _("apply a view transform based on personal or camera maker look,\n"
       "for corrective purposes, to prepare images for display"),
     _("corrective"),
     _("linear, RGB, display-referred"),
     _("non-linear, RGB"),
     _("non-linear, RGB, display-referred"));
}

int default_group()
{
  return IOP_GROUP_BASIC | IOP_GROUP_TECHNICAL;
}

int flags()
{
  return IOP_FLAGS_SUPPORTS_BLENDING | IOP_FLAGS_ALLOW_TILING;
}

dt_iop_colorspace_type_t default_colorspace(dt_iop_module_t *self,
                                            dt_dev_pixelpipe_t *pipe,
                                            dt_dev_pixelpipe_iop_t *piece)
{
  return IOP_CS_RGB;
}

static void set_presets(dt_iop_module_so_t *self,
                        const basecurvergb_preset_t *presets,
                        const int count,
                        const gboolean camera)
{
  dt_develop_blend_params_t default_blendop_params;
  dt_develop_blend_init_blend_parameters(&default_blendop_params, DEVELOP_BLEND_CS_RGB_DISPLAY);

  // transform presets above to db entries
  for(int k = 0; k < count; k++)
  {
    dt_iop_basecurvergb_params_t tmp = presets[k].params;
    gchar *prefixed_name = camera ? g_strdup(presets[k].name) : g_strdup_printf(BUILTIN_PREFIX "%s", presets[k].name);
    // add the preset.
    dt_gui_presets_add_with_blendop(prefixed_name, self->op, self->version(),
                                    &tmp, sizeof(dt_iop_basecurvergb_params_t),
                                    &default_blendop_params, 1);
    // and restrict it to model, maker, iso, and raw images
    dt_gui_presets_update_mml(prefixed_name, self->op, self->version(),
                              presets[k].maker, presets[k].model, "");
    dt_gui_presets_update_iso(prefixed_name, self->op, self->version(),
                              presets[k].iso_min, presets[k].iso_max);
    dt_gui_presets_update_format(prefixed_name, self->op, self->version(), FOR_RAW);
    // make it auto-apply for matching images:
    dt_gui_presets_update_autoapply(prefixed_name, self->op, self->version(), FALSE);
    // hide all non-matching presets in case the model string is set.
    // When force_autoapply was given always filter (as these are per-camera presets)
    dt_gui_presets_update_filter(prefixed_name,
                                 self->op, self->version(), camera || presets[k].filter);
    g_free(prefixed_name);
  }
}

static gboolean _match(const char *value, const char *pattern)
{
  char *pat = g_strdup(pattern);

  // the pattern is for SQL, replace '%' by '*' and '_' by '.'

  int k=0;
  while(pat[k])
  {
    if(pat[k] == '%')
      pat[k] = '*';
    else if(pat[k] == '_')
      pat[k] = '.';
    k++;
  }

  gboolean res = g_regex_match_simple(pat, value, G_REGEX_CASELESS, G_REGEX_MATCH_ANCHORED);

  g_free(pat);

  return res;
}

static gboolean _check_camera(dt_iop_basecurvergb_params_t *d,
                              const char *e_maker,
                              const char *e_model,
                              const char *c_maker,
                              const char *c_model,
                              const basecurvergb_preset_t *presets,
                              const int count)
{
  // in reverse order as the more specific maker/models is after
  // the more generic and we want to match the more specific.
  for(int k = count - 1; k > 0; k--)
  {
    if((_match(e_maker, presets[k].maker)
        && _match(e_model, presets[k].model))
       || (_match(c_maker, presets[k].maker)
           && _match(c_model, presets[k].model)))
    {
      *d = presets[k].params;
      return TRUE;
    }
  }

  return FALSE;
}

void reload_defaults(dt_iop_module_t *self)
{
  dt_iop_basecurvergb_params_t *const d = self->default_params;

  if(self->multi_priority == 0)
  {
    const dt_image_t *const image = &(self->dev->image_storage);

    self->default_enabled = FALSE;

    gboolean FOUND = FALSE;

    // first check for camera specific basecure if needed

    const gboolean autoapply_percamera =
      dt_conf_get_bool("plugins/darkroom/basecurve/auto_apply_percamera_presets");

    if(autoapply_percamera)
    {
      FOUND = _check_camera(d,
                            image->exif_maker, image->exif_model,
                            image->camera_maker, image->camera_alias,
                            basecurvergb_camera_presets, basecurvergb_camera_presets_cnt);
    }

    if(!FOUND)
    {
      // then check for default base curve for the camera

      FOUND = _check_camera(d,
                            image->exif_maker, image->exif_model,
                            image->camera_maker, image->camera_alias,
                            basecurvergb_presets, basecurvergb_presets_cnt);
    }
  }
  else
  {
    // set to neutral (cubic-spline) for all other instances
    *d = basecurvergb_presets[0].params;
    d->preserve_hue = 1.0f;
    d->preserve_highlight_saturation = 0.0f;
    d->source_white = 0.0f;
  }
}

void init_presets(dt_iop_module_so_t *self)
{
  // sql begin
  dt_database_start_transaction(darktable.db);

  set_presets(self, basecurvergb_presets, basecurvergb_presets_cnt, FALSE);
  set_presets(self, basecurvergb_camera_presets, basecurvergb_camera_presets_cnt, TRUE);

  // sql commit
  dt_database_release_transaction(darktable.db);

  // auto-applied display-referred default
  self->pref_based_presets = TRUE;

  const gboolean is_display_referred = dt_is_display_referred();

  if(is_display_referred)
  {
    dt_gui_presets_add_generic
      (_("display-referred default"), self->op, self->version(),
       NULL, 0,
       1, DEVELOP_BLEND_CS_RGB_DISPLAY);

    dt_gui_presets_update_format(BUILTIN_PRESET("display-referred default"), self->op,
                                 self->version(), FOR_RAW);

    dt_gui_presets_update_autoapply(BUILTIN_PRESET("display-referred default"),
                                    self->op, self->version(), TRUE);
  }
}


#if 0
void tiling_callback(dt_iop_module_t *self,
                     dt_dev_pixelpipe_iop_t *piece,
                     const dt_iop_roi_t *roi_in,
                     const dt_iop_roi_t *roi_out,
                     dt_develop_tiling_t *tiling)
{
    tiling->factor = 2.0f;                   // in + out
    tiling->maxbuf = 1.0f;
    tiling->overhead = 0;
    tiling->xalign = 1;
    tiling->yalign = 1;
    tiling->overlap = 0;
}
#endif

void process(dt_iop_module_t *self,
             dt_dev_pixelpipe_iop_t *piece,
             const void *const ivoid,
             void *const ovoid,
             const dt_iop_roi_t *const roi_in,
             const dt_iop_roi_t *const roi_out)
{

  const dt_iop_order_iccprofile_info_t *const work_profile
      = dt_ioppr_get_pipe_current_profile_info(self, piece->pipe);
  if(work_profile == NULL) return; // cannot continue without a working profile


  const float *const in = (const float *)ivoid;
  float *const out = (float *)ovoid;
  //const int ch = piece->colors; <-- it appears someone was trying to make this handle monochrome data,
  //however the for loops only handled RGBA - FIXME, determine what possible data formats and channel
  //configurations we might encounter here and handle those too
  dt_iop_basecurvergb_data_t *const d = piece->data;
  // TODO: will be used for color preservation
  // const dt_iop_order_iccprofile_info_t *const work_profile = dt_ioppr_get_iop_work_profile_info(piece->module, piece->module->dev->iop);

  const int wd = roi_in->width, ht = roi_in->height;
  const float factor_source_white = powf(2.0f, d->source_white);
  const float preserve_hue = d->preserve_hue;


  // get matrix for working profile to XYZ_D65 conversion
  // TODO: I assume that working profile matrix is conversion to XYZ_D50, so that adaption to D65 is needed. Correct?
  dt_colormatrix_t XYZ_D65_to_working_profile = { { 0.0f } };
  dt_colormatrix_t working_profile_to_XYZ_D65 = { { 0.0f } };
  dt_colormatrix_mul(working_profile_to_XYZ_D65, XYZ_D50_to_D65_CAT16, work_profile->matrix_in);
  dt_colormatrix_mul(XYZ_D65_to_working_profile, work_profile->matrix_out, XYZ_D65_to_D50_CAT16);

  // make the transposed matrices
  dt_colormatrix_t XYZ_D65_to_working_profile_transposed;
  dt_colormatrix_t working_profile_to_XYZ_D65_transposed;
  dt_colormatrix_transpose(XYZ_D65_to_working_profile_transposed, XYZ_D65_to_working_profile);
  dt_colormatrix_transpose(working_profile_to_XYZ_D65_transposed, working_profile_to_XYZ_D65);


  // TODO: scale to user-defined max value and clip
  const size_t npixels = (size_t)wd * ht;
  DT_OMP_FOR()
  for(size_t k = 0; k < 4*npixels; k += 4)
  {
    for(int i = 0; i < 3; i++)
    {
      const float in_multiplied = in[k+i] * factor_source_white;
      // use base curve for values < 1, else use extrapolation.
      if(in_multiplied < 1.0f)
        out[k+i] = fmaxf(d->table[CLAMP((int)(in_multiplied * 0x10000ul), 0, 0xffff)], 0.f);
      else
        out[k+i] = fmaxf(dt_iop_eval_exp(d->unbounded_coeffs, in_multiplied), 0.f);
    }

    // TODO: apply preserve_hue here
    dt_aligned_pixel_t RGB_in;
    dt_aligned_pixel_t RGB_out;
    dt_aligned_pixel_t RGB_out_before_hue_preservation;
    copy_pixel(RGB_in, in + k);
    copy_pixel(RGB_out, out + k);
    copy_pixel(RGB_out_before_hue_preservation, out + k);

    dt_aligned_pixel_t XYZ_D65_in;
    dt_aligned_pixel_t XYZ_D65_out;
    dt_apply_transposed_color_matrix(RGB_in, working_profile_to_XYZ_D65_transposed, XYZ_D65_in);
    dt_apply_transposed_color_matrix(RGB_out, working_profile_to_XYZ_D65_transposed, XYZ_D65_out);

    // to jab (todo: later we will use oklab for preserve_hue)
    dt_aligned_pixel_t jab_in;
    dt_aligned_pixel_t jab_out;
    dt_XYZ_D65_to_oklab(XYZ_D65_in, jab_in);
    dt_XYZ_D65_to_oklab(XYZ_D65_out, jab_out);

    dt_aligned_pixel_t jch_in;
    dt_aligned_pixel_t jch_out;
    // dt_JzAzBz_2_JzCzhz is general JAB to JCH conversion and can be used
    dt_JzAzBz_2_JzCzhz(jab_in, jch_in);
    dt_JzAzBz_2_JzCzhz(jab_out, jch_out);

    // insert hue from in to out
    jch_out[2] = preserve_hue * jch_in[2] + (1.0f - preserve_hue) * jch_out[2];

    // convert back to working profile
    dt_JzCzhz_2_JzAzBz(jch_out, jab_out);
    dt_oklab_to_XYZ_D65(jab_out, XYZ_D65_out);

    // convert back to RGB
    dt_apply_transposed_color_matrix(XYZ_D65_out, XYZ_D65_to_working_profile_transposed, RGB_out);

    // TODO: recover lightness after highlight saturation preservation
    // by comparing lightness before and after saturation preservation
    // dt_aligned_pixel_t HSL_out_before_hue_preservation;
    // dt_aligned_pixel_t HSL_out_after_hue_preservation;
    // dt_RGB_2_HSV(RGB_out_before_hue_preservation, HSL_out_before_hue_preservation);
    // dt_RGB_2_HSV(RGB_out, HSL_out_after_hue_preservation);
    // const float L_before = HSL_out_before_hue_preservation[2];
    // // const float L_after = HSL_out_after_hue_preservation[1];
    // HSL_out_after_hue_preservation[2] = L_before * 3;
    // dt_HSV_2_RGB(HSL_out_after_hue_preservation, RGB_out);


    // do saturation preservation
    const float min = fminf(RGB_out[0], fminf(RGB_out[1], RGB_out[2]));
    const float max = fmaxf(RGB_out[0], fmaxf(RGB_out[1], RGB_out[2]));
    const float delta = max - min;

    const float L = (min + max) / 2.0f;
    float C;

    if(fabsf(max) > 1e-6f && fabsf(delta) > 1e-6f)
    {
      C = delta;
    }
    else
    {
      C = 0.0f;
    }
    const float factor_resaturation = sqrtf(L * C) * d->preserve_highlight_saturation;

    dt_aligned_pixel_t HSV_out;
    dt_aligned_pixel_t HSV_in;
    dt_RGB_2_HSV(RGB_in, HSV_in);
    dt_RGB_2_HSV(RGB_out, HSV_out);
    HSV_out[1] = HSV_in[1] * factor_resaturation + (1.0f - factor_resaturation) * HSV_out[1];
    dt_HSV_2_RGB(HSV_out, RGB_out);

    out[k+0] = RGB_out[0];
    out[k+1] = RGB_out[1];
    out[k+2] = RGB_out[2];    
    // out[k+0] = factor_resaturation;
    // out[k+1] = factor_resaturation;
    // out[k+2] = factor_resaturation;
    out[k+3] = in[k+3];


    // TODO: see colorbalancergb on how to get work profile
    // like   dt_colormatrix_mul(output_matrix, XYZ_D50_to_D65_CAT16, work_profile->matrix_in); // output_matrix used as temp buffer
    // and then to get the color matrix etc. for converting from work profile to XYZ_D65, which is needed for oklab conversion

    // Convert in and out to oklab

    // copy hue from in to out

    // Convert back to RGB


    // TODO: apply preserve_hue here
  }
}

void commit_params(dt_iop_module_t *self,
                   dt_iop_params_t *p1,
                   dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_basecurvergb_data_t *d = piece->data;
  dt_iop_basecurvergb_params_t *p = (dt_iop_basecurvergb_params_t *)p1;

  d->preserve_hue = p->preserve_hue;
  d->preserve_highlight_saturation = p->preserve_highlight_saturation;
  d->source_white = p->source_white;

  const int ch = 0;
  // take care of possible change of curve type or number of nodes (not yet implemented in UI)
  if(d->basecurvergb_type != p->basecurvergb_type[ch] || d->basecurvergb_nodes != p->basecurvergb_nodes[ch])
  {
    if(d->curve) // catch initial init_pipe case
      dt_draw_curve_destroy(d->curve);
    d->curve = dt_draw_curve_new(0.0, 1.0, p->basecurvergb_type[ch]);
    d->basecurvergb_nodes = p->basecurvergb_nodes[ch];
    d->basecurvergb_type = p->basecurvergb_type[ch];
    for(int k = 0; k < p->basecurvergb_nodes[ch]; k++)
    {
      // printf("p->basecurve[%i][%i].x = %f;\n", ch, k, p->basecurve[ch][k].x);
      // printf("p->basecurve[%i][%i].y = %f;\n", ch, k, p->basecurve[ch][k].y);
      (void)dt_draw_curve_add_point(d->curve, p->basecurve[ch][k].x, p->basecurve[ch][k].y);
    }
  }
  else
  {
    for(int k = 0; k < p->basecurvergb_nodes[ch]; k++)
      dt_draw_curve_set_point(d->curve, k, p->basecurve[ch][k].x, p->basecurve[ch][k].y);
  }
  dt_draw_curve_calc_values(d->curve, 0.0f, 1.0f, 0x10000, NULL, d->table);

  // now the extrapolation stuff:
  const float xm = p->basecurve[0][p->basecurvergb_nodes[0] - 1].x;
  const float x[4] = { 0.7f * xm, 0.8f * xm, 0.9f * xm, 1.0f * xm };
  const float y[4] = { d->table[CLAMP((int)(x[0] * 0x10000ul), 0, 0xffff)],
                       d->table[CLAMP((int)(x[1] * 0x10000ul), 0, 0xffff)],
                       d->table[CLAMP((int)(x[2] * 0x10000ul), 0, 0xffff)],
                       d->table[CLAMP((int)(x[3] * 0x10000ul), 0, 0xffff)] };
  dt_iop_estimate_exp(x, y, 4, d->unbounded_coeffs);
}

void init_pipe(dt_iop_module_t *self,
               dt_dev_pixelpipe_t *pipe,
               dt_dev_pixelpipe_iop_t *piece)
{
  // create part of the pixelpipe
  piece->data = calloc(1, sizeof(dt_iop_basecurvergb_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(dt_iop_module_t *self,
                  dt_dev_pixelpipe_t *pipe,
                  dt_dev_pixelpipe_iop_t *piece)
{
  // clean up everything again.
  dt_iop_basecurvergb_data_t *d = piece->data;
  dt_draw_curve_destroy(d->curve);
  free(piece->data);
  piece->data = NULL;
}

void gui_update(dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;

  // gui curve is read directly from params during expose event.
  gtk_widget_queue_draw(GTK_WIDGET(g->area));
}

static float eval_grey(float x)
{
  // "log base" is a combined scaling and offset change so that x->[0,1], with
  // the left side of the histogram expanded (slider->right) or not (slider left, linear)
  return x;
}

void init(dt_iop_module_t *self)
{
  dt_iop_default_init(self);
  dt_iop_basecurvergb_params_t *d = self->default_params;
  d->basecurve[0][1].x = d->basecurve[0][1].y = 1.0;
  d->basecurvergb_nodes[0] = 2;
}

void init_global(dt_iop_module_so_t *self)
{
  dt_iop_basecurvergb_global_data_t *gd = malloc(sizeof(dt_iop_basecurvergb_global_data_t));
  self->data = gd;
}

void cleanup_global(dt_iop_module_so_t *self)
{
  free(self->data);
  self->data = NULL;
}

static gboolean dt_iop_basecurvergb_leave_notify(GtkWidget *widget,
                                              GdkEventCrossing *event,
                                              dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  if(!(event->state & GDK_BUTTON1_MASK))
    g->selected = -1;
  gtk_widget_queue_draw(widget);
  return FALSE;
}

/**
 * Applies log scaling to input value x based on the base parameter.
 * 
 * @param x Input value in range [0, 1]
 * @param base Scaling base:
 *             > 0: Spreads shadows (left side) - larger values compress more towards highlights
 *             < 0: Spreads highlights (right side) - larger absolute values compress more towards shadows  
 *             = 0: Linear mapping (no transformation)
 * @return Transformed value in range [0, 1]
 */
static float to_log(const float x, const float base)
{
  if(base > 0.0f)
    return logf(x * base + 1.0f) / logf(base + 1.0f);
  else if(base < 0.0f)
  {
    // For negative base values, spread highlights by applying log transform to (1-x) 
    // and then mirroring the result around 0.5
    const float abs_base = -base;
    const float flipped_x = 1.0f - x;
    const float log_result = logf(flipped_x * abs_base + 1.0f) / logf(abs_base + 1.0f);
    return 1.0f - log_result;
  }
  else
    return x;
}

/**
 * Applies the inverse of log scaling - converts from log-scaled space back to linear.
 * This is the mathematical inverse of to_log() function.
 * 
 * @param x Input value in log-scaled space [0, 1]  
 * @param base Same base parameter used in corresponding to_log() call
 * @return Linear value in range [0, 1]
 */
static float to_lin(const float x, const float base)
{
  if(base > 0.0f)
    return (powf(base + 1.0f, x) - 1.0f) / base;
  else if(base < 0.0f)
  {
    // Inverse transformation for negative base values
    // Mirror x around 0.5, apply inverse log transform, then mirror back
    const float abs_base = -base;
    const float flipped_x = 1.0f - x;
    const float linear_result = (powf(abs_base + 1.0f, flipped_x) - 1.0f) / abs_base;
    return 1.0f - linear_result;
  }
  else
    return x;
}

static gboolean dt_iop_basecurvergb_draw(GtkWidget *widget, cairo_t *crf, dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  dt_iop_basecurvergb_params_t *p = self->params;

  int nodes = p->basecurvergb_nodes[0];
  dt_iop_basecurvergb_node_t *basecurve = p->basecurve[0];
  if(g->minmax_curve_type != p->basecurvergb_type[0] || g->minmax_curve_nodes != p->basecurvergb_nodes[0])
  {
    dt_draw_curve_destroy(g->minmax_curve);
    g->minmax_curve = dt_draw_curve_new(0.0, 1.0, p->basecurvergb_type[0]);
    g->minmax_curve_nodes = p->basecurvergb_nodes[0];
    g->minmax_curve_type = p->basecurvergb_type[0];
    for(int k = 0; k < p->basecurvergb_nodes[0]; k++)
      (void)dt_draw_curve_add_point(g->minmax_curve, p->basecurve[0][k].x, p->basecurve[0][k].y);
  }
  else
  {
    for(int k = 0; k < p->basecurvergb_nodes[0]; k++)
      dt_draw_curve_set_point(g->minmax_curve, k, p->basecurve[0][k].x, p->basecurve[0][k].y);
  }
  dt_draw_curve_t *minmax_curve = g->minmax_curve;
  dt_draw_curve_calc_values(minmax_curve, 0.0, 1.0, DT_IOP_TONECURVE_RES, g->draw_xs, g->draw_ys);

  float unbounded_coeffs[3];
  const float xm = basecurve[nodes - 1].x;
  {
    const float x[4] = { 0.7f * xm, 0.8f * xm, 0.9f * xm, 1.0f * xm };
    const float y[4] = { g->draw_ys[CLAMP((int)(x[0] * DT_IOP_TONECURVE_RES), 0, DT_IOP_TONECURVE_RES - 1)],
                         g->draw_ys[CLAMP((int)(x[1] * DT_IOP_TONECURVE_RES), 0, DT_IOP_TONECURVE_RES - 1)],
                         g->draw_ys[CLAMP((int)(x[2] * DT_IOP_TONECURVE_RES), 0, DT_IOP_TONECURVE_RES - 1)],
                         g->draw_ys[CLAMP((int)(x[3] * DT_IOP_TONECURVE_RES), 0, DT_IOP_TONECURVE_RES - 1)] };
    dt_iop_estimate_exp(x, y, 4, unbounded_coeffs);
  }

  const int inset = DT_GUI_CURVE_EDITOR_INSET;
  GtkAllocation allocation;
  gtk_widget_get_allocation(widget, &allocation);
  int width = allocation.width, height = allocation.height;
  cairo_surface_t *cst = dt_cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create(cst);
  // clear bg
  cairo_set_source_rgb(cr, .2, .2, .2);
  cairo_paint(cr);

  cairo_translate(cr, inset, inset);
  width -= 2 * inset;
  height -= 2 * inset;

  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.0));
  cairo_set_source_rgb(cr, .1, .1, .1);
  cairo_rectangle(cr, 0, 0, width, height);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, .3, .3, .3);
  cairo_rectangle(cr, 0, 0, width, height);
  cairo_fill(cr);

  cairo_translate(cr, 0, height);
  if(g->selected >= 0)
  {
    char text[30];
    // draw information about current selected node
    PangoLayout *layout;
    PangoRectangle ink;
    PangoFontDescription *desc = pango_font_description_copy_static(darktable.bauhaus->pango_font_desc);
    pango_font_description_set_weight(desc, PANGO_WEIGHT_BOLD);
    pango_font_description_set_absolute_size(desc, PANGO_SCALE);
    layout = pango_cairo_create_layout(cr);
    pango_layout_set_font_description(layout, desc);

    const float x_node_value = basecurve[g->selected].x * 100;
    const float y_node_value = basecurve[g->selected].y * 100;
    const float d_node_value = y_node_value - x_node_value;
    // scale conservatively to 100% of width:
    snprintf(text, sizeof(text), "100.00 / 100.00 ( +100.00)");
    pango_layout_set_text(layout, text, -1);
    pango_layout_get_pixel_extents(layout, &ink, NULL);
    pango_font_description_set_absolute_size(desc, (double)width / ink.width * PANGO_SCALE);
    pango_layout_set_font_description(layout, desc);

    snprintf(text, sizeof(text), "%.2f / %.2f ( %+.2f)", x_node_value, y_node_value, d_node_value);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    pango_layout_set_text(layout, text, -1);
    pango_layout_get_pixel_extents(layout, &ink, NULL);
    cairo_move_to(cr, 0.98f * width - ink.width - ink.x, -0.02 * height - ink.height - ink.y);
    pango_cairo_show_layout(cr, layout);
    cairo_stroke(cr);
    pango_font_description_free(desc);
    g_object_unref(layout);
  }
  cairo_scale(cr, 1.0f, -1.0f);

  // draw grid
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(.4));
  cairo_set_source_rgb(cr, .1, .1, .1);
  if(g->loglogscale != 0.0f)
  {
    // Custom grid drawing that matches our to_log function behavior
    const int num = 4;
    for(int k = 1; k < num; k++)
    {
      const float grid_pos = k / (float)num;
      const float x = to_log(grid_pos, g->loglogscale);
      const float y = to_log(grid_pos, g->loglogscale);
      
      // Vertical lines
      cairo_move_to(cr, x * width, 0);
      cairo_line_to(cr, x * width, height);
      cairo_stroke(cr);
      
      // Horizontal lines  
      cairo_move_to(cr, 0, y * height);
      cairo_line_to(cr, width, y * height);
      cairo_stroke(cr);
    }
  }
  else
    dt_draw_grid(cr, 4, 0, 0, width, height);

  // draw nodes positions
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.));
  cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
  for(int k = 0; k < nodes; k++)
  {
    const float x = to_log(basecurve[k].x, g->loglogscale), y = to_log(basecurve[k].y, g->loglogscale);
    cairo_arc(cr, x * width, y * height, DT_PIXEL_APPLY_DPI(3), 0, 2. * M_PI);
    cairo_stroke(cr);
  }

  // draw selected cursor
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.));

  if(g->selected >= 0)
  {
    cairo_set_source_rgb(cr, .9, .9, .9);
    const float x = to_log(basecurve[g->selected].x, g->loglogscale),
                y = to_log(basecurve[g->selected].y, g->loglogscale);
    cairo_arc(cr, x * width, y * height, DT_PIXEL_APPLY_DPI(4), 0, 2. * M_PI);
    cairo_stroke(cr);
  }

  // draw curve
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(2.));
  cairo_set_source_rgb(cr, .9, .9, .9);
  // cairo_set_line_cap  (cr, CAIRO_LINE_CAP_SQUARE);
  cairo_move_to(cr, 0, height * to_log(g->draw_ys[0], g->loglogscale));
  for(int k = 1; k < DT_IOP_TONECURVE_RES; k++)
  {
    const float xx = k / (DT_IOP_TONECURVE_RES - 1.0f);
    if(xx > xm)
    {
      const float yy = dt_iop_eval_exp(unbounded_coeffs, xx);
      const float x = to_log(xx, g->loglogscale), y = to_log(yy, g->loglogscale);
      cairo_line_to(cr, x * width, height * y);
    }
    else
    {
      const float yy = g->draw_ys[k];
      const float x = to_log(xx, g->loglogscale), y = to_log(yy, g->loglogscale);
      cairo_line_to(cr, x * width, height * y);
    }
  }
  cairo_stroke(cr);

  cairo_destroy(cr);
  cairo_set_source_surface(crf, cst, 0, 0);
  cairo_paint(crf);
  cairo_surface_destroy(cst);
  return TRUE;
}

static inline int _add_node(dt_iop_basecurvergb_node_t *basecurve, int *nodes, float x, float y)
{
  int selected = -1;
  if(basecurve[0].x > x)
    selected = 0;
  else
  {
    for(int k = 1; k < *nodes; k++)
    {
      if(basecurve[k].x > x)
      {
        selected = k;
        break;
      }
    }
  }
  if(selected == -1) selected = *nodes;
  for(int i = *nodes; i > selected; i--)
  {
    basecurve[i].x = basecurve[i - 1].x;
    basecurve[i].y = basecurve[i - 1].y;
  }
  // found a new point
  basecurve[selected].x = x;
  basecurve[selected].y = y;
  (*nodes)++;
  return selected;
}

static void dt_iop_basecurvergb_sanity_check(dt_iop_module_t *self, GtkWidget *widget)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  dt_iop_basecurvergb_params_t *p = self->params;

  int ch = 0;
  int nodes = p->basecurvergb_nodes[ch];
  dt_iop_basecurvergb_node_t *basecurve = p->basecurve[ch];

  if(nodes <= 2) return;

  const float mx = basecurve[g->selected].x;

  // delete vertex if order has changed
  // for all points, x coordinate of point must be strictly larger than
  // the x coordinate of the previous point
  if((g->selected > 0 && (basecurve[g->selected - 1].x >= mx))
     || (g->selected < nodes - 1 && (basecurve[g->selected + 1].x <= mx)))
  {
    for(int k = g->selected; k < nodes - 1; k++)
    {
      basecurve[k].x = basecurve[k + 1].x;
      basecurve[k].y = basecurve[k + 1].y;
    }
    g->selected = -2; // avoid re-insertion of that point immediately after this
    p->basecurvergb_nodes[ch]--;
  }
}

static gboolean _move_point_internal(dt_iop_module_t *self,
                                     GtkWidget *widget,
                                     float dx,
                                     float dy,
                                     const guint state);

static gboolean dt_iop_basecurvergb_motion_notify(GtkWidget *widget,
                                               GdkEventMotion *event,
                                               dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  dt_iop_basecurvergb_params_t *p = self->params;
  int ch = 0;
  int nodes = p->basecurvergb_nodes[ch];
  dt_iop_basecurvergb_node_t *basecurve = p->basecurve[ch];

  GtkAllocation allocation;
  gtk_widget_get_allocation(widget, &allocation);
  const int inset = DT_GUI_CURVE_EDITOR_INSET;
  int height = allocation.height - 2 * inset, width = allocation.width - 2 * inset;
  const double old_m_x = g->mouse_x;
  const double old_m_y = g->mouse_y;
  g->mouse_x = event->x - inset;
  g->mouse_y = event->y - inset;

  const float mx = CLAMP(g->mouse_x, 0, width) / (float)width;
  const float my = 1.0f - CLAMP(g->mouse_y, 0, height) / (float)height;
  const float linx = to_lin(mx, g->loglogscale), liny = to_lin(my, g->loglogscale);

  if(event->state & GDK_BUTTON1_MASK)
  {
    // got a vertex selected:
    if(g->selected >= 0)
    {
      // this is used to translate mause position in loglogscale to make this behavior unified with linear scale.
      const float translate_mouse_x = old_m_x / width - to_log(basecurve[g->selected].x, g->loglogscale);
      const float translate_mouse_y = 1 - old_m_y / height - to_log(basecurve[g->selected].y, g->loglogscale);
      // dx & dy are in linear coordinates
      const float dx = to_lin(g->mouse_x / width - translate_mouse_x, g->loglogscale)
                       - to_lin(old_m_x / width - translate_mouse_x, g->loglogscale);
      const float dy = to_lin(1 - g->mouse_y / height - translate_mouse_y, g->loglogscale)
                       - to_lin(1 - old_m_y / height - translate_mouse_y, g->loglogscale);

      return _move_point_internal(self, widget, dx, dy, event->state);
    }
    else if(nodes < MAXNODES && g->selected >= -1)
    {
      // no vertex was close, create a new one!
      g->selected = _add_node(basecurve, &p->basecurvergb_nodes[ch], linx, liny);
      dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
    }
  }
  else
  {
    // minimum area around the node to select it:
    float min = .04f;
    min *= min; // comparing against square
    int nearest = -1;
    for(int k = 0; k < nodes; k++)
    {
      float dist
          = (my - to_log(basecurve[k].y, g->loglogscale)) * (my - to_log(basecurve[k].y, g->loglogscale))
            + (mx - to_log(basecurve[k].x, g->loglogscale)) * (mx - to_log(basecurve[k].x, g->loglogscale));
      if(dist < min)
      {
        min = dist;
        nearest = k;
      }
    }
    g->selected = nearest;
  }
  if(g->selected >= 0) gtk_widget_grab_focus(widget);
  gtk_widget_queue_draw(widget);
  return TRUE;
}

static gboolean dt_iop_basecurvergb_button_press(GtkWidget *widget,
                                              GdkEventButton *event,
                                              dt_iop_module_t *self)
{
  dt_iop_basecurvergb_params_t *p = self->params;
  const dt_iop_basecurvergb_params_t *const d = self->default_params;
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;

  int ch = 0;
  int nodes = p->basecurvergb_nodes[ch];
  dt_iop_basecurvergb_node_t *basecurve = p->basecurve[ch];

  if(event->button == GDK_BUTTON_PRIMARY)
  {
    if(event->type == GDK_BUTTON_PRESS && dt_modifier_is(event->state, GDK_CONTROL_MASK)
      && nodes < MAXNODES && g->selected == -1)
    {
      // if we are not on a node -> add a new node at the current x of the pointer and y of the curve at that x
      const int inset = DT_GUI_CURVE_EDITOR_INSET;
      GtkAllocation allocation;
      gtk_widget_get_allocation(widget, &allocation);
      int width = allocation.width - 2 * inset;
      g->mouse_x = event->x - inset;
      g->mouse_y = event->y - inset;

      const float mx = CLAMP(g->mouse_x, 0, width) / (float)width;
      const float linx = to_lin(mx, g->loglogscale);

      // don't add a node too close to others in x direction, it can crash dt
      int selected = -1;
      if(basecurve[0].x > linx)
        selected = 0;
      else
      {
        for(int k = 1; k < nodes; k++)
        {
          if(basecurve[k].x > linx)
          {
            selected = k;
            break;
          }
        }
      }
      if(selected == -1) selected = nodes;
      // > 0 -> check distance to left neighbour
      // < nodes -> check distance to right neighbour
      if(!((selected > 0 && linx - basecurve[selected - 1].x <= 0.025) ||
           (selected < nodes && basecurve[selected].x - linx <= 0.025)))
      {
        // evaluate the curve at the current x position
        const float y = dt_draw_curve_calc_value(g->minmax_curve, linx);

        if(y >= 0.0 && y <= 1.0) // never add something outside the viewport, you couldn't change it afterwards
        {
          // create a new node
          selected = _add_node(basecurve, &p->basecurvergb_nodes[ch], linx, y);

          // maybe set the new one as being selected
          float min = .04f;
          min *= min; // comparing against square
          for(int k = 0; k < nodes; k++)
          {
            float other_y = to_log(basecurve[k].y, g->loglogscale);
            float dist = (y - other_y) * (y - other_y);
            if(dist < min) g->selected = selected;
          }

          dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
          gtk_widget_queue_draw(GTK_WIDGET(g->area));
        }
      }
      return TRUE;
    }
    else if(event->type == GDK_2BUTTON_PRESS)
    {
      // reset current curve
      p->basecurvergb_nodes[ch] = d->basecurvergb_nodes[ch];
      p->basecurvergb_type[ch] = d->basecurvergb_type[ch];
      for(int k = 0; k < d->basecurvergb_nodes[ch]; k++)
      {
        p->basecurve[ch][k].x = d->basecurve[ch][k].x;
        p->basecurve[ch][k].y = d->basecurve[ch][k].y;
      }
      g->selected = -2; // avoid motion notify re-inserting immediately.
      dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
      gtk_widget_queue_draw(GTK_WIDGET(g->area));
      return TRUE;
    }
  }
  else if(event->button == GDK_BUTTON_SECONDARY && g->selected >= 0)
  {
    if(g->selected == 0 || g->selected == nodes - 1)
    {
      float reset_value = g->selected == 0 ? 0 : 1;
      basecurve[g->selected].y = basecurve[g->selected].x = reset_value;
      gtk_widget_queue_draw(GTK_WIDGET(g->area));
      dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
      return TRUE;
    }

    for(int k = g->selected; k < nodes - 1; k++)
    {
      basecurve[k].x = basecurve[k + 1].x;
      basecurve[k].y = basecurve[k + 1].y;
    }
    basecurve[nodes - 1].x = basecurve[nodes - 1].y = 0;
    g->selected = -2; // avoid re-insertion of that point immediately after this
    p->basecurvergb_nodes[ch]--;
    gtk_widget_queue_draw(GTK_WIDGET(g->area));
    dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
    return TRUE;
  }
  return FALSE;
}

static gboolean _move_point_internal(dt_iop_module_t *self,
                                     GtkWidget *widget,
                                     float dx,
                                     float dy,
                                     const guint state)
{
  dt_iop_basecurvergb_params_t *p = self->params;
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;

  const int ch = 0;
  dt_iop_basecurvergb_node_t *basecurve = p->basecurve[ch];

  const float multiplier = dt_accel_get_speed_multiplier(widget, state);
  dx *= multiplier;
  dy *= multiplier;

  basecurve[g->selected].x = CLAMP(basecurve[g->selected].x + dx, 0.0f, 1.0f);
  basecurve[g->selected].y = CLAMP(basecurve[g->selected].y + dy, 0.0f, 1.0f);

  dt_iop_basecurvergb_sanity_check(self, widget);

  gtk_widget_queue_draw(widget);
  dt_dev_add_history_item_target(darktable.develop, self, TRUE, widget);
  return TRUE;
}

#define basecurvergb_DEFAULT_STEP (0.001f)

static gboolean _scrolled(GtkWidget *widget, GdkEventScroll *event, dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;

  if(dt_gui_ignore_scroll(event)) return FALSE;

  if(g->selected < 0) return TRUE;

  gdouble delta_y;
  if(dt_gui_get_scroll_delta(event, &delta_y))
  {
    delta_y *= -basecurvergb_DEFAULT_STEP;
    return _move_point_internal(self, widget, 0.0, delta_y, event->state);
  }

  return TRUE;
}

static gboolean dt_iop_basecurvergb_key_press(GtkWidget *widget,
                                           GdkEventKey *event,
                                           dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;

  if(g->selected < 0) return TRUE;

  int handled = 0;
  float dx = 0.0f, dy = 0.0f;
  if(event->keyval == GDK_KEY_Up || event->keyval == GDK_KEY_KP_Up)
  {
    handled = 1;
    dy = basecurvergb_DEFAULT_STEP;
  }
  else if(event->keyval == GDK_KEY_Down || event->keyval == GDK_KEY_KP_Down)
  {
    handled = 1;
    dy = -basecurvergb_DEFAULT_STEP;
  }
  else if(event->keyval == GDK_KEY_Right || event->keyval == GDK_KEY_KP_Right)
  {
    handled = 1;
    dx = basecurvergb_DEFAULT_STEP;
  }
  else if(event->keyval == GDK_KEY_Left || event->keyval == GDK_KEY_KP_Left)
  {
    handled = 1;
    dx = -basecurvergb_DEFAULT_STEP;
  }

  if(!handled) return FALSE;

  return _move_point_internal(self, widget, dx, dy, event->state);
}

#undef basecurvergb_DEFAULT_STEP

void gui_changed(dt_iop_module_t *self, GtkWidget *w, void *previous)
{
}

static void logbase_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  g->loglogscale = eval_grey(dt_bauhaus_slider_get(g->logbase));
  gtk_widget_queue_draw(GTK_WIDGET(g->area));
}

void gui_init(dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = IOP_GUI_ALLOC(basecurvergb);
  const dt_iop_basecurvergb_params_t *const p = self->default_params;

  g->minmax_curve = dt_draw_curve_new(0.0, 1.0, p->basecurvergb_type[0]);
  g->minmax_curve_type = p->basecurvergb_type[0];
  g->minmax_curve_nodes = p->basecurvergb_nodes[0];
  for(int k = 0; k < p->basecurvergb_nodes[0]; k++)
    (void)dt_draw_curve_add_point(g->minmax_curve, p->basecurve[0][k].x, p->basecurve[0][k].y);
  g->mouse_x = g->mouse_y = -1.0;
  g->selected = -1;
  g->loglogscale = 0;

  g->area = GTK_DRAWING_AREA(dtgtk_drawing_area_new_with_height(0));
  gtk_widget_set_tooltip_text(GTK_WIDGET(g->area), _("abscissa: input, ordinate: output. works on RGB channels"));
  g_object_set_data(G_OBJECT(g->area), "iop-instance", self);
  dt_action_define_iop(self, NULL, N_("curve"), GTK_WIDGET(g->area), NULL);

  self->widget = dt_gui_vbox(g->area);

  // initially set to 1 (consistency with previous versions), but double-click resets to 0
  // to get a quick way to reach 0 with the mouse.
  g->preserve_hue = dt_bauhaus_slider_from_params(self, "preserve_hue");
  dt_bauhaus_slider_set_default(g->preserve_hue, 0.0f);
  dt_bauhaus_slider_set_digits(g->preserve_hue, 3);
  gtk_widget_set_tooltip_text(g->preserve_hue, _("The strength of hue preservation after application of base curve"));
  gtk_widget_set_no_show_all(g->preserve_hue, TRUE);
  gtk_widget_set_visible(g->preserve_hue, CL_TRUE);

  g->preserve_highlight_saturation = dt_bauhaus_slider_from_params(self, "preserve_highlight_saturation");
  dt_bauhaus_slider_set_default(g->preserve_highlight_saturation, 0.0f);
  dt_bauhaus_slider_set_digits(g->preserve_highlight_saturation, 3);
  gtk_widget_set_tooltip_text(g->preserve_highlight_saturation, _("The strength of hue preservation after application of base curve"));
  gtk_widget_set_no_show_all(g->preserve_highlight_saturation, TRUE);
  gtk_widget_set_visible(g->preserve_highlight_saturation, CL_TRUE);

  g->source_white = dt_bauhaus_slider_from_params(self, "source_white");
  dt_bauhaus_slider_set_default(g->source_white, 1.0f);
  dt_bauhaus_slider_set_digits(g->source_white, 3);
  gtk_widget_set_tooltip_text(g->source_white, _("Number of ev stops the source white lies over / below 1.0"));
  gtk_widget_set_no_show_all(g->source_white, TRUE);
  gtk_widget_set_visible(g->source_white, CL_TRUE);

  g->logbase = dt_bauhaus_slider_new_with_range(self, -40.0f, 40.0f, 0, 0.0f, 2);
  dt_bauhaus_widget_set_label(g->logbase, NULL, N_("scale for graph"));
  g_signal_connect(G_OBJECT(g->logbase), "value-changed", G_CALLBACK(logbase_callback), self);
  dt_gui_box_add(self->widget, g->logbase);

  gtk_widget_add_events(GTK_WIDGET(g->area), GDK_POINTER_MOTION_MASK | darktable.gui->scroll_mask
                                           | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK
                                           | GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK);
  gtk_widget_set_can_focus(GTK_WIDGET(g->area), TRUE);
  g_signal_connect(G_OBJECT(g->area), "draw", G_CALLBACK(dt_iop_basecurvergb_draw), self);
  g_signal_connect(G_OBJECT(g->area), "button-press-event", G_CALLBACK(dt_iop_basecurvergb_button_press), self);
  g_signal_connect(G_OBJECT(g->area), "motion-notify-event", G_CALLBACK(dt_iop_basecurvergb_motion_notify), self);
  g_signal_connect(G_OBJECT(g->area), "leave-notify-event", G_CALLBACK(dt_iop_basecurvergb_leave_notify), self);
  g_signal_connect(G_OBJECT(g->area), "scroll-event", G_CALLBACK(_scrolled), self);
  g_signal_connect(G_OBJECT(g->area), "key-press-event", G_CALLBACK(dt_iop_basecurvergb_key_press), self);
}

void gui_cleanup(dt_iop_module_t *self)
{
  dt_iop_basecurvergb_gui_data_t *g = self->gui_data;
  dt_draw_curve_destroy(g->minmax_curve);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
