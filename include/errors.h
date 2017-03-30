/*
canal is C program to calculate electrical conductivities from
molecular simulation data.
Copyright 2017 Frank Uhlig (uhlig.frank@gmail.com)

This file is part of canal.

canal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

canal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with canal.  If not, see <http://www.gnu.org/licenses/>.
*/


typedef enum errors_e
{
  ERROR=1000,
  FILE_NOT_FOUND,
  MISSING_INPUT_PARAM,
  NOT_IMPLEMENTED,
  NONSENSICAL,
  INCOMPLETE_INPUT,
  NOT_ENOUGH_DATA,
  FATAL,
  EXPCODE,
  UNASSIGNED_ERROR,
  IO_ERROR,

} errors_t;


typedef enum warnings_e
{
  WARNING=2000,
  MEMORY_WASTE,
  YOU_KNOW_WHAT,
  NO_EFFECT,
  EXPWARNING,
  EXPFEATURE,
  UNASSIGNED_WARNING,

} warnings_t;

void print_error(int errno, char * detail, char * file, int line);
void print_warning(int errno, char * detail);
void print_error_header();
void print_error_footer();
void print_warning_header();
void print_warning_footer();
