@echo off
cls
color 1B
mode con:cols=98
mode con:lines=50
if not exist PhyML-3.1_win32.exe goto error
if exist PhyML-3.1_win32.exe then goto launch
echo on


:launch
PhyML-3.1_win32
goto end

:error
echo Error - can't find `PhyML-3.1_win32.exe'
pause
goto end


:end
echo Execution finished
pause
