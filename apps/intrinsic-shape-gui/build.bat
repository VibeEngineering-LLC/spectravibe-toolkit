@echo off
rem Build exe (one-folder) into dist\IntrinsicShapeGUI\. Python: env ISGUI_PYTHON (venv with requirements.txt), else python from PATH.
chcp 65001 >nul
cd /d "%~dp0"
if "%ISGUI_PYTHON%"=="" set ISGUI_PYTHON=python
rem version 0.2.0-beta.<sha>: short sha of HEAD goes into the build (src\isgui\_build_sha.txt, not in git)
git rev-parse --short HEAD > src\isgui\_build_sha.txt 2>nul || echo dev> src\isgui\_build_sha.txt
"%ISGUI_PYTHON%" -m PyInstaller IntrinsicShapeGUI.spec --noconfirm --clean
