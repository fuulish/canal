#include "errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_error(int errno, char * detail, char * file, int line)
{

  print_error_header();
  
  switch ( errno )
  {
  case FILE_NOT_FOUND:
    printf("File '%s' could not be opened\n", detail);
    break;
  case MISSING_INPUT_PARAM:
    printf("Task cannot be completed without specifying input parameter: '%s'\n", detail);
    break;
  case NOT_IMPLEMENTED:
    printf("%s not implemented!\n", detail);
    break;
  case NONSENSICAL:
    printf("%s is non-sensical!\n", detail);
    break;
  case INCOMPLETE_INPUT:
    printf("%s\n", detail);
    break;
  case NOT_ENOUGH_DATA:
    printf("Not enough data in '%s'\n", detail);
    break;
  case FATAL:
    printf("Fatal error: %s\n", detail);
    break;
  case EXPCODE:
    printf("%s is only available in experimental mode. Re-configure with --enable-experimental.\n", detail);
    break;
  case IO_ERROR:
    printf("%s IO problem detected.\n", detail);
    break;
  case UNASSIGNED_ERROR:
  default:
    printf("An error occurred that has not yet been assigned to a category!\n");
    break;
  }

  printf( "Stopping in file %s at line %d\n", file, line );
  
  print_error_footer();
  
  exit ( errno );
}

void print_warning(int warno, char * detail)
{

  print_warning_header();
  
  switch ( warno )
  {
  case MEMORY_WASTE:
    printf("We are wasting memory, because %s\n", detail);
    break;
  case YOU_KNOW_WHAT:
    printf("%s\n", detail);
    break;
  case NO_EFFECT:
    printf("%s has no effect.\n", detail);
    break;
  case EXPWARNING:
    printf("Compiled with experimental features (just a reminder ;-).\n");
    break;
  case EXPFEATURE:
    printf("Using experimental feature: %s\n", detail);
    break;
  case UNASSIGNED_WARNING:
  default:
    printf("You are doing something possibly dangerous, but no warning has been assigned yet. Sorry!\n");
    break;
    
  }
  
  print_warning_footer();

}

void print_error_header()
{
  printf("\n!!!!!!!!! E R R O R !!!!!!!!!\n\n");
}

void print_error_footer()
{
  printf("\n!!!!!!!!! E R R O R !!!!!!!!!\n\n");
}

void print_warning_header()
{
  printf("\n!!!!!!!!! W A R N I N G !!!!!!!!!\n\n");
}

void print_warning_footer()
{
  printf("\n!!!!!!!!! W A R N I N G !!!!!!!!!\n\n");
}
