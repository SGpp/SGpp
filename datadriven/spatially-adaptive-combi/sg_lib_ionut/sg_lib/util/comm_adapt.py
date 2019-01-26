import numpy as np
from math import trunc

class CommAdapt(object):

	def __truncate(self, number, digits):

		stepper = pow(10.0, digits)
		res     = trunc(stepper * number) / stepper

		return res

	def read_real_part(self, file_name):

		real_part = 0.
		with open(file_name, "r") as file:
		    lines           = file.readlines()
		    data_line       = lines[-1]

		    data_split = data_line.split()
		    if len(data_split) == 2:
		        real_part_str, imag_part_str    = data_split
		        real_part                       = self.__truncate(float(real_part_str), 6)
		    elif len(data_split) == 3:
		        ky, real_part_str, imag_part_str        = data_split
		        real_part                               = self.__truncate(float(real_part_str), 6)

		file.close()

		return real_part


	def read_imag_part(self, file_name):

		imag_part = 0.
		with open(file_name, "r") as file:
			lines 		= file.readlines()
			data_line 	= lines[-1]

			data_split = data_line.split()
			if len(data_split) == 2:
				real_part_str, imag_part_str 	= data_split
				imag_part 						= float(imag_part_str)
			elif len(data_split) == 3:
				ky, real_part_str, imag_part_str 	= data_split
				imag_part 							= float(imag_part_str)

		file.close()

		return imag_part

	def write_sg_data(self, file_name, sg_data):

		with open(file_name, "w") as file:
			for sg_info in sg_data:
				line = ''

				sg_points 	= sg_info[0]
				sim_no 		= sg_info[1]

				for sg_point in sg_points:
					line += str(sg_point) + ' '

				line += str(sim_no) + '\n'
				
				file.write(line)

		file.close

	def write_is_terminated(self, file_name, is_not_terminated):
		
		is_not_terminated = str(is_not_terminated).lower()

		with open(file_name, "w") as file:
			file.write(is_not_terminated + '\n')

		file.close()