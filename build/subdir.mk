# BSD 2-Clause License
#
# Copyright (c) [2019] [Valeria Barra]
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

################################################################################
# For any file added to the project, this file needs to be updated
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../QR.f90 \
../domain_space.f90 \
../domain_time.f90 \
../get_der_cap_x.f90 \
../get_der_hhh_x.f90 \
../get_der_quadterm.f90 \
../get_der_quadvdw_x.f90 \
../get_der_vdw_x.f90 \
../get_matrix_x.f90 \
../get_matrix_x_c.f90 \
../get_matrix_x_hxave.f90 \
../get_matrix_x_quadterm.f90 \
../get_matrix_x_quadvdw.f90 \
../get_matrix_x_vdw.f90 \
../hhh_cap.f90 \
../init.f90 \
../main.f90 \
../mobility.f90 \
../new_time.f90 \
../newton.f90 \
../nonlin_cap.f90 \
../nonlin_func.f90 \
../nonlin_h.f90 \
../nonlin_quadterm.f90 \
../nonlin_quadvdw.f90 \
../nonlin_vdw.f90 \
../nr_ban.f90 \
../nrtype.f90 \
../output.f90 \
../par.f90 \
../paras.f90 \
../penta.f90 \
../term_c.f90 \
../term_coeff_one.f90 \
../term_coeff_three.f90 \
../term_quad.f90 \
../term_quadvdw.f90 \
../term_vdw.f90 \
../vanderwaals.f90 

OBJS += \
./QR.o \
./domain_space.o \
./domain_time.o \
./get_der_cap_x.o \
./get_der_hhh_x.o \
./get_der_quadterm.o \
./get_der_quadvdw_x.o \
./get_der_vdw_x.o \
./get_matrix_x2.o \
./get_matrix_x_c.o \
./get_matrix_x_hxave.o \
./get_matrix_x_quadterm.o \
./get_matrix_x_quadvdw.o \
./get_matrix_x_vdw.o \
./hhh_cap.o \
./init.o \
./main.o \
./mobility.o \
./new_time.o \
./newton.o \
./nonlin_cap.o \
./nonlin_func.o \
./nonlin_h.o \
./nonlin_quadterm.o \
./nonlin_quadvdw.o \
./nonlin_vdw.o \
./nr_ban.o \
./nrtype.o \
./output.o \
./par.o \
./paras.o \
./penta.o \
./term_c.o \
./term_coeff_one.o \
./term_coeff_three.o \
./term_quad.o \
./term_quadvdw.o \
./term_vdw.o \
./vanderwaals.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

QR.o: ../QR.f90 domain_space.o nrtype.o paras.o

domain_space.o: ../domain_space.f90 nrtype.o

domain_time.o: ../domain_time.f90 nrtype.o

get_der_cap_x.o: ../get_der_cap_x.f90 domain_space.o mobility.o nrtype.o paras.o

get_der_hhh_x.o: ../get_der_hhh_x.f90 domain_space.o nrtype.o paras.o

get_der_quadterm.o: ../get_der_quadterm.f90 domain_space.o mobility.o nrtype.o paras.o

get_der_quadvdw_x.o: ../get_der_quadvdw_x.f90 domain_space.o nrtype.o paras.o vanderwaals.o

get_der_vdw_x.o: ../get_der_vdw_x.f90 domain_space.o nrtype.o paras.o vanderwaals.o

get_matrix_x.o: ../get_matrix_x.f90 domain_space.o domain_time.o nrtype.o paras.o

get_matrix_x_c.o: ../get_matrix_x_c.f90 domain_space.o nrtype.o paras.o

get_matrix_x_hxave.o: ../get_matrix_x_hxave.f90 domain_space.o nrtype.o paras.o

get_matrix_x_quadterm.o: ../get_matrix_x_quadterm.f90 domain_space.o nrtype.o paras.o

get_matrix_x_quadvdw.o: ../get_matrix_x_quadvdw.f90 domain_space.o nrtype.o paras.o

get_matrix_x_vdw.o: ../get_matrix_x_vdw.f90 domain_space.o nrtype.o paras.o

hhh_cap.o: ../hhh_cap.f90 domain_space.o nrtype.o

init.o: ../init.f90 domain_space.o domain_time.o nrtype.o paras.o 

main.o: ../main.f90 domain_space.o domain_time.o nrtype.o paras.o

mobility.o: ../mobility.f90 nrtype.o

new_time.o: ../new_time.f90 domain_space.o domain_time.o mobility.o nr_ban.o nrtype.o paras.o vanderwaals.o

newton.o: ../newton.f90 domain_space.o domain_time.o mobility.o nr_ban.o nrtype.o paras.o vanderwaals.o

nonlin_cap.o: ../nonlin_cap.f90 domain_space.o mobility.o nrtype.o

nonlin_func.o: ../nonlin_func.f90 domain_space.o mobility.o nrtype.o paras.o vanderwaals.o

nonlin_h.o: ../nonlin_h.f90 domain_space.o nrtype.o

nonlin_quadterm.o: ../nonlin_quadterm.f90 domain_space.o mobility.o nrtype.o

nonlin_quadvdw.o: ../nonlin_quadvdw.f90 domain_space.o nrtype.o vanderwaals.o

nonlin_vdw.o: ../nonlin_vdw.f90 domain_space.o mobility.o nrtype.o vanderwaals.o

nr_ban.o: ../nr_ban.f90 nrtype.o

nrtype.o: ../nrtype.f90

output.o: ../output.f90 domain_space.o domain_time.o nrtype.o paras.o

par.o: ../par.f90 domain_space.o paras.o

paras.o: ../paras.f90 nrtype.o

penta.o: ../penta.f90 nr_ban.o nrtype.o

term_c.o: ../term_c.f90 domain_space.o nrtype.o paras.o

term_cgravitycos.o: ../term_cgravitycos.f90 domain_space.o nrtype.o paras.o

term_coeff_three.o: ../term_coeff_three.f90 domain_space.o nrtype.o paras.o

term_quad.o: ../term_quad.f90 domain_space.o nrtype.o paras.o

term_quadvdw.o: ../term_quadvdw.f90 domain_space.o nrtype.o paras.o

term_vdw.o: ../term_vdw.f90 domain_space.o nrtype.o paras.o

vanderwaals.o: ../vanderwaals.f90 nrtype.o paras.o


