
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
