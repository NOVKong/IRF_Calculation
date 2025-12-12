@ECHO OFF
setlocal EnableDelayedExpansion
::
IF {"%~1"}=={""} ( set "F=L_ULS_Sd.txt"
) else set "F=%~1"
IF {"%~2"}=={""} ( set /a N=10
) else set /a N=%~2
echo 1st Arg:%1; 2nd Arg:%2 

set "D=!F:~,-4!"
if not exist "!D!\r" md "!D!\r"
if exist !F! copy/y !F! !D!\
:: Not to store visuallation...
if exist ..\Wisting_L.json copy ..\Wisting_L.json !D!\.
:: commands = sre_textcmd(LC_FILE,startat=START_Line)
::TxtRep startat=START_Line startat=0
if exist ..\runSRE0.py TxtRep LC_FILE "'!F!'" ..\runSRE0.py !D!\runSRE.py

if exist "!D!\allrun_!N!.bat" del /q "!D!\allrun_!N!.bat" 2>nul
set nLC=1
for /f "usebackq delims=" %%A in ("%F%") do (
 echo ===^>!D!\r\!nLC!.bat
    echo %%A>"!D!\r\!nLC!.bat"
    set /a mod=!nLC! %% !N!
    if !mod! == 0 (
        echo set timestamp=!N!:%%date%% %%time%%>>"!D!\allrun_!N!.bat"
        echo title run !D! - !N!%%timestamp%% >>"!D!\allrun_!N!.bat"
        echo cmd /c r\!nLC!.bat>>"!D!\allrun_!N!.bat"
    ) else  echo start "!nLC!" cmd /c r\!nLC!.bat>>"!D!\allrun_!N!.bat"
    set /a nLC+=1
)
goto :EOF

SimaPp /eD:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\B_ULS_Sd\LC10x3h\rLC2\node_1\key_sima_elmfor.txt /lelm.txt /bI /d0
SimaPp /nD:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\B_ULS_Sd\LC10x3h\rLC2\node_1\key_sima_noddis.txt /lnod.txt /bI /d0
SimaPp /rD:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\B_ULS_Sd\LC10x3h\rLC2\node_1\results.txt /lres.txt /bI /d0

for %D in (B_ULS_Sd.txt,B_ALS1F_Sd100yr.txt,B_ALS2F_Sd10yr.txt,B_ULS_Sd10kyr.txt,B_ALS1F_Sd1yr.txt,B_ALS2F_Sd1yr.txt,B_ALS3F_Sd1yr.txt) do batchGenB %D 20
for %D in (L_ULS_Sd.txt,L_ALS1F_Sd100yr.txt,L_ALS2F_Sd10yr.txt,L_ULS_Sd10kyr.txt,L_ALS1F_Sd1yr.txt,L_ALS2F_Sd1yr.txt,L_ALS3F_Sd1yr.txt) do batchGenL %D 20

for %d in (B_ULS_Sd,B_ALS1F_Sd100yr,B_ALS2F_Sd10yr,B_ULS_Sd10kyr,B_ALS1F_Sd1yr,B_ALS2F_Sd1yr,B_ALS3F_Sd1yr) do (title %d&..\mv_simaRuns.bat %d 0 DST=E:\WORK\M2114_Wisting\SIMA\CStudy)
for %d in (L_ULS_Sd,L_ALS1F_Sd100yr,L_ALS2F_Sd10yr,L_ULS_Sd10kyr,L_ALS1F_Sd1yr,L_ALS2F_Sd1yr,L_ALS3F_Sd1yr) do (title %d&..\mv_simaRuns.bat %d 0 DST=E:\WORK\M2114_Wisting\SIMA\CStudy)

D:\WorkDir\ProjSpace\M2114_Wisting\Sima>python find_missing_lcs.py 
start "B_ULS_Sd10kyr" cmd /c "timeout 20800&Python runSre.py"

for %d in (B_ULS_Sd,B_ALS1F_Sd100yr,B_ALS2F_Sd10yr,B_ULS_Sd10kyr,B_ALS1F_Sd1yr,B_ALS2F_Sd1yr,B_ALS3F_Sd1yr) do (title %d&pushd %d&python runSRE.py&popd)
for %d in (L_ULS_Sd,L_ALS1F_Sd100yr,L_ALS2F_Sd10yr,L_ULS_Sd10kyr,L_ALS1F_Sd1yr,L_ALS2F_Sd1yr,L_ALS3F_Sd1yr) do (title %d&pushd %d&python runSRE.py&popd)

REM Sreening
D:\WorkDir\ProjSpace\M2114_Wisting\Sima>mv_simaRuns.bat CStudy\OpenWater_Loaded\L_ULS_Screen
D:\WorkDir\ProjSpace\M2114_Wisting\Sima>ScanLC_Error.bat CStudy\OpenWater_Loaded\L_ULS_Screen\r
for %d in (B_ULS_Screen) do (python ..\run_post_parallel.py %d --hour=-1)
for %d in (B_ULS_Screen) do (python ..\run_simappos.py %d -i=-1 -m MTF)


for %d in (B_ULS_Sd,B_ALS1F_Sd100yr,B_ALS2F_Sd10yr,B_ULS_Sd10kyr,B_ALS1F_Sd1yr,B_ALS2F_Sd1yr,B_ALS3F_Sd1yr) do (python ..\run_post_parallel.py %d --hour=-1)
for %d in (B_ULS_Sd,B_ALS1F_Sd100yr,B_ALS2F_Sd10yr,B_ULS_Sd10kyr,B_ALS1F_Sd1yr,B_ALS2F_Sd1yr,B_ALS3F_Sd1yr) do (python ..\run_simappos.py %d -3=3 -i=-1 -m Obs)

for %d in (L_ULS_Sd,L_ALS1F_Sd100yr,L_ALS2F_Sd10yr,L_ULS_Sd10kyr,L_ALS1F_Sd1yr,L_ALS2F_Sd1yr,L_ALS3F_Sd1yr) do (python ..\run_post_parallel.py %d --hour=-1)
for %d in (L_ULS_Sd,L_ALS1F_Sd100yr,L_ALS2F_Sd10yr,L_ULS_Sd10kyr,L_ALS1F_Sd1yr,L_ALS2F_Sd1yr,L_ALS3F_Sd1yr) do (python ..\run_simappos.py %d -3=3 -i=-1 -m Obs)

=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd100yr.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd10yr.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ULS_Sd10kyr.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS1F_Sd1yr.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS2F_Sd1yr.xlsx]SIMAPlt'!$AK$5:$AK$16)
=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$T$5:$T$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$AG$5:$AG$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$AH$5:$AH$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$AI$5:$AI$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$AJ$5:$AJ$16)	=MAX('D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Ballast\[B_ALS3F_Sd1yr.xlsx]SIMAPlt'!$AK$5:$AK$16)


python ..\run_post_parallel.py D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Loaded\L_ALS3F_Sd1yr\rLC46\node_1 --hour=-1
ModelTest2Xls /p"D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Loaded\L_ALS3F_Sd1yr\rLC46\node_1" "/lLALS3F_1yr\rLC46_000303002_070806_14" /x"L_ALS3F_Sd1yr.xlsx" /t"SIMA SIMAPlt" /sX /MObs /M3hour /MBL=19748 /c1000.0 /n12

python ..\run_post_parallel.py D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Loaded\L_ALS1F_Sd1yr\rLC10\node_1 --hour=-1
ModelTest2Xls /p"D:\WorkDir\ProjSpace\M2114_Wisting\Sima\SRE_LowMLBE_Loaded\L_ALS1F_Sd1yr\rLC10\node_1" "/lLALS1F_1yr\rLC10_000303002_07_11" /x"L_ALS1F_Sd1yr.xlsx" /t"SIMA SIMAPlt" /sX /MObs /M3hour /MBL=19748 /c1000.0 /n12



robocopy "C:\Users\abc\AppData\Local" "E:\WORK\AppData\Local" /mir /xj
rmdir "C:\Users\YourName\AppData\Local" /s /q
mklink /j "C:\Users\YourName\AppData\Local" "E:\WORK\AppData\Local"