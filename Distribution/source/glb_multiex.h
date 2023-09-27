/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef GLB_MULTIEX_H
#define GLB_MULTIEX_H 1

#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif

#include "glb_types.h"
#include "glb_error.h"
/* Options for glbInitExpFromParent */
#define GLB_COPY_RULES      1
#define GLB_DONT_COPY_RULES 0
#define GLB_FAST_RATES 1
#define GLB_SLOW_RATES 0


extern int glb_ignore_invalid_chi2;
extern int glb_num_of_rules;
extern int glb_current_experiment;
glb_nuisance *glb_alloc_nuisance();
int glb_copy_nuisance(glb_nuisance *dest, glb_nuisance *src);
int glb_free_nuisance(glb_nuisance *n);

glb_exp glbAllocExp();
int glbDefaultExp(glb_exp ins);
void glbInitExp(glb_exp ins);
void glbInitExpFromParent(struct glb_experiment *exp, struct glb_experiment *p);
int glbSetRatesInExperiment(int exp, int which_rates, int fast_rates);
void glbResetExp(struct glb_experiment *in);
void glbFreeExp(struct glb_experiment *in);
void glbExpAddChild(struct glb_experiment *parent, struct glb_experiment *child);
void glbExpRemoveChild(struct glb_experiment *parent, struct glb_experiment *child);
int glbPrintExpByPointer(struct glb_experiment *exp);

int glbRateTemplate(struct glb_experiment *e, int which_rates);

void glb_set_profile_scaling(double scale,int i);

/* Declarations for functions defined in the AEDL parser */
glb_naming *glb_copy_names_from_parser(glb_naming *head);
glb_naming *glb_copy_names(glb_naming *in, glb_naming *head);
void glb_free_names(glb_naming *stale);

#endif /* GLB_MULTIEX_H */
