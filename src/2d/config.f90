module config
#include "../config.h"

    ! Name of package
    character(*), parameter :: package = PACKAGE

    ! Define to the address where bug reports for this package should be sent.
    character(*), parameter :: package_bugreport = PACKAGE_BUGREPORT

    ! Define to the full name of this package.
    character(*), parameter :: package_name = PACKAGE_NAME

    ! Define to the full name and version of this package.
    character(*), parameter :: package_string = PACKAGE_STRING

    ! Define to the one symbol short name of this package.
    character(*), parameter :: package_tarname = PACKAGE_TARNAME

    ! Define to the home page for this package.
    character(*), parameter :: package_url = PACKAGE_URL

    ! Define to the version of this package.
    character(*), parameter :: package_version = PACKAGE_VERSION

    ! Version number of package
    character(*), parameter :: version = VERSION

end module config
