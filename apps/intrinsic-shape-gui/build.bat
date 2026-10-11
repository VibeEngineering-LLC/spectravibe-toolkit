@echo off
rem Build exe (one-folder) into dist\IntrinsicShapeGUI\. Python: env ISGUI_PYTHON (venv with requirements.txt), else python from PATH.
chcp 65001 >nul
cd /d "%~dp0"
if "%ISGUI_PYTHON%"=="" set ISGUI_PYTHON=python
"%ISGUI_PYTHON%" -m PyInstaller IntrinsicShapeGUI.spec --noconfirm --clean
