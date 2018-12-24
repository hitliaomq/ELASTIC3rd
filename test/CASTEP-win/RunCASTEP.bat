@echo off
rem ---------------------------------------------------------------------------
rem
rem Script for stand-alone CASTEP execution.
rem
rem This program and all subroutines, data, and files used by it
rem are protected by Copyright and hence may not be used, copied,
rem modified, transmitted, inspected, or executed by any means including
rem the use of electronic data processing equipment, xerography, or
rem any other methods without the express written permission of the
rem copyright holder.
rem
rem Copyright (c) 2015, Dassault Systemes, All Rights Reserved
rem
rem ***************************************************************************
SETLOCAL

set MS_INSTALL_ROOT=C:\Program Files (x86)\BIOVIA\Materials Studio 17.2

set server=Castep
call "%MS_INSTALL_ROOT%"\share\bin\runMSserver.bat %server% %*
if ERRORLEVEL 6 goto usage

goto end

:usage
type  "%MS_INSTALL_ROOT%"\etc\CASTEP\bin\RunCASTEP.Readme

:end

ENDLOCAL
