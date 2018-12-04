@echo off
if "%1" == "dll" (
setlocal
#SET _REQUIRED_INCLUDE_PATHS=/ID:\Apps\AMD\acml-6.1.0\ifort64\include
#SET _REQUIRED_LIBS=D:\Apps\AMD\acml-6.1.0\ifort64\lib\libacml_dll.lib
SET _REQUIRED_INCLUDE_PATHS=/ID:\Apps\fftw\include
SET _REQUIRED_LIBS=D:\Apps\fftw\lib\fftw3.lib
ECHO Required include paths: %_REQUIRED_INCLUDE_PATHS%
ECHO Required libraries: %_REQUIRED_LIBS%
cl /D_USRDLL /D_WINDLL /LD /DBUILDING_GABOR_DLL %_REQUIRED_INCLUDE_PATHS% /Fo:build\gabor.o include\gabor.c %_REQUIRED_LIBS% /link /DLL /IMPLIB:.\lib\libgabor.lib /out:.\lib\libgabor.dll
endlocal
)
if "%1" == "python" (
setlocal EnableDelayedExpansion
SET REQUIRED_INCLUDE_PATH=/ID:\Apps\fftw\include /ID:\Apps\Anaconda3\include /ID:\Apps\Anaconda3\Lib\site-packages\numpy\core\include
SET REQUIRED_LIBS=D:\Apps\fftw\lib\fftw3.lib D:\Apps\Anaconda3\libs\python36.lib
cl /LD  /DBUILDING_PYHTON_MODULE /DNDEBUG !REQUIRED_INCLUDE_PATH! /Fo:build\gabor.o include\gabor.c !REQUIRED_LIBS! /link /DLL /IMPLIB:.\lib\libgabor.lib /out:.\lib\gabor.pyd
endlocal
if "%2" == "install" (
echo Copy module to: %PYTHONPATH%\
copy .\lib\gabor.pyd %PYTHONPATH%\
)
)
