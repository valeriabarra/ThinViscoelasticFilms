# BSD 2-Clause (MIT) License
#
# Copyright (c) [2019] [ThinViscoelasticFilms]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
../get_matrix_x2.f90 \
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
../nonlin_func2.f90 \
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
./nonlin_func2.o \
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

get_matrix_x2.o: ../get_matrix_x2.f90 domain_space.o domain_time.o nrtype.o paras.o

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

nonlin_func2.o: ../nonlin_func2.f90 domain_space.o mobility.o nrtype.o paras.o vanderwaals.o

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


