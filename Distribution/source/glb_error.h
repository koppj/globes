/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2004,  The GLoBES Team
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




#ifndef GLB_ERROR_H
#define GLB_ERROR_H 1

#include <stdlib.h>
#include <stdio.h>

extern const char *glb_prog_name;
extern void glb_prog_name_init (const char *argv0);

extern int glbSetVerbosityLevel(int level);


extern void glb_warning      (const char *message);
extern void glb_error        (const char *message);
extern void glb_fatal        (const char *message);

extern void *glb_malloc (size_t size);
extern void *glb_realloc (void *ptr, size_t size);
extern void glb_free(void *ptr);
extern FILE *glb_fopen(const char *filename, const char *mode);
extern int glb_fclose(FILE *stream);

#endif /* !GLB_ERROR_H */
