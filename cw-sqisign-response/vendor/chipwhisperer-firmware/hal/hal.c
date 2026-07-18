/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2015 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hal.h"

__attribute__((weak)) void led_ok(unsigned int status)
{
}

__attribute__((weak)) void led_error(unsigned int status)
{
}

#ifdef __GNUC__
#if ((__GNUC__ > 11) || \
     ((__GNUC__ == 11) && (__GNUC_MINOR__ >= 3)))
__attribute__((weak)) void _close() {}
__attribute__((weak)) void _fstat() {}
__attribute__((weak)) void _getpid() {}
__attribute__((weak)) void _isatty() {}
__attribute__((weak)) void _kill() {}
__attribute__((weak)) void _lseek() {}
__attribute__((weak)) void _read() {}
__attribute__((weak)) void _write() {}
#endif
#endif