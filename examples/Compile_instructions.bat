if "%1" == "fourier" cl /Fo:.\build\test_fouriertransform.obj /ID:\Apps\fftw\include .\src\test_fouriertransform.c D:\Apps\fftw\lib\fftw3.lib /link /OUT:.\bin\test_fouriertransform.exe
if "%1" == "gabor"  cl /Fo:.\build\test_gaborfilter.obj /I..\include /IC:\fftw\include .\src\test_gaborfilter.c C:\fftw\lib\fftw3.lib ..\lib\libgabor.lib /link /OUT:.\bin\test_gaborfilter.exe
