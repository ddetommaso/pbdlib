### PbDLib C++ Library ###

http://programming-by-demonstration.org


Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh, Leonel Rozo, Tohid Alizadeh

    PbDLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PbDLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PbDLib.  If not, see <http://www.gnu.org/licenses/>.



# Contributors:
    Davide De Tommaso ----> dtmdvd[at]gmail[dot]com
    Milad Malekzadeh  ----> milad[dot]malekzadeh[at]gmail[dot]com
    Leonel Rozo ----------> ing[dot]leonelrozo[at]gmail[dot]com
    Tohid Alizadeh -------> tohid[dot]alizadeh[at]gmail[dot]com



# 1. Compiling instructions


    ### 1.1 Dependences (you can find in deps directory)
	(1) Armadillo C++ versions >3.9
	(2) cmake

    ### 1.2 Building

	$ cd directory_of_the_pbdlib
	$ mkdir build
	$ cd build
	$ cmake ..
	$ make

    ### 1.3 Testing
	$ cd directory_of_the_pbdlib
	$ cd build/test
	$ ./test_datapoints
