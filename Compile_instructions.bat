setlocal
#SET _REQUIRED_INCLUDE_PATHS=/ID:\Apps\AMD\acml-6.1.0\ifort64\include
#SET _REQUIRED_LIBS=D:\Apps\AMD\acml-6.1.0\ifort64\lib\libacml_dll.lib
SET _REQUIRED_INCLUDE_PATHS=/IC:\fftw\include
SET _REQUIRED_LIBS=C:\fftw\lib\fftw3.lib
ECHO Required include paths: %_REQUIRED_INCLUDE_PATHS%
ECHO Required libraries: %_REQUIRED_LIBS%
cl /D_USRDLL /D_WINDLL /LD /DBUILDING_GABOR_DLL %_REQUIRED_INCLUDE_PATHS% /Fo:build\gabor.o include\gabor.c %_REQUIRED_LIBS% /link /DLL /IMPLIB:.\lib\libgabor.lib /out:.\lib\libgabor.dll
endlocal
if "%1" == "examples" (
copy .\lib\libgabor.dll .\examples
)
