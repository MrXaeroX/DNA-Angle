@echo off
:: ------------------------------------------------------------------------------
:: Copyright (C) 2020 Alexander V. Popov.
::
:: This source code is free software; you can redistribute it and/or
:: modify it under the terms of the GNU General Public License as
:: published by the Free Software Foundation; either version 2 of
:: the License, or (at your option) any later version.
::
:: This source code is distributed in the hope that it will be
:: useful, but WITHOUT ANY WARRANTY; without even the implied
:: warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
:: See the GNU General Public License for more details.
::
:: You should have received a copy of the GNU General Public License
:: along with this program; if not, write to the Free Software
:: Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
:: ------------------------------------------------------------------------------
setlocal
set project_dir=%~dp0
set pause_on_exit=0
echo %cmdcmdline% | find /i "%~0" >nul
if not errorlevel 1 set pause_on_exit=1
set PYTHONDONTWRITEBYTECODE=1

:: Defer control.
python "%project_dir%\dna_angle.py" %*
if _%pause_on_exit%_==_1_ pause
