### PbDLib C++ Library ###

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


# 1. What's PbDLib ?

    PbDLib is an open source C++ library for using the Programming-by-Demonstration machine learning tools in your C++ code. Most of the tools provided by PbDLib are implementations of the algorithms described in the book "Robot Programming by Demonstration: A Probabilistic Approach" (Sylvain Calinon, 2009) and other related scientific publications.

    You can download the source code at http://github.com/ddetommaso/pbdlib

    For more information about Programming-by-Demonstration please visit http://programming-by-demonstration.org


# 2. Contributors:

    PbDLib has been developed as part of the research activities inside the Learning and Interaction Group, Advanced Robotics Dept at the Istituto Italiano di Tecnologia, Genova.
    In the following you can find the list of the contributors

    Davide De Tommaso ----> dtmdvd[at]gmail[dot]com
    Milad Malekzadeh  ----> milad[dot]malekzadeh[at]gmail[dot]com
    Leonel Rozo ----------> ing[dot]leonelrozo[at]gmail[dot]com
    Tohid Alizadeh -------> tohid[dot]alizadeh[at]gmail[dot]com

    For more information about the research of the group please visit
    http://www.iit.it/en/advr-labs/learning-and-interaction.html


# 3. Compiling instructions


    ### 1.1 Dependences (you can find in deps directory)
	(1) Armadillo C++ versions >3.9
	(2) CMake
	(3) Doxygen

    ### 1.2 Building
	$ mkdir build
	$ cd build
	$ cmake ..
	$ make

    ### 1.3 Testing
	$ cd build/test
	$ ./test_datapoints

    ### 1.4 Documentation generation
	$ cd doc
	$ doxygen Doxyfile


# 4. Some useful citations

    Did you find PbDLib useful for your research?

    Please consider to acknowledge the authors in any academic publications that have made use of this code or part of it. In the following you can find some code related BibTex references:


    @book{Calinon09book,
	author="S. Calinon",
	year="2009",
	title="Robot Programming by Demonstration: A Probabilistic Approach",
	publisher="EPFL/CRC Press"
	note="{EPFL} {P}ress {ISBN} 978-2-940222-31-5, {CRC} {P}ress {ISBN} 978-1-4398-0867-2"
    }


    @article{Calinon07SMC,
  	title="On Learning, Representing and Generalizing a Task in a Humanoid Robot",
  	author="S. Calinon and F. Guenter and A. Billard",
  	journal="IEEE Transactions on Systems, Man and Cybernetics, Part B. Special issue on robot learning by observation, demonstration and imitation",
  	year="2007",
  	volume="37",
  	number="2",
  	pages="286--298"
    }


    @inproceedings{Calinon14ICRA,
	author="Calinon, S. and Bruno, D. and Caldwell, D. G.",
	title="A task-parameterized probabilistic model with minimal intervention control",
	booktitle="Proc. {IEEE} Intl Conf. on Robotics and Automation ({ICRA})",
	year="2014",
	month="May-June",
	address="Hong Kong, China",
	pages="3339--3344"
    }


    @inproceedings{Calinon12Hum,
	author="Calinon, S. and Li, Z. and Alizadeh, T. and Tsagarakis, N. G. and Caldwell, D. G.",
	title="Statistical dynamical systems for skills acquisition in humanoids",
	booktitle="Proc. {IEEE} Intl Conf. on Humanoid Robots ({H}umanoids)",
	year="2012",
	address="Osaka, Japan",
	pages="323--329"
    }

# 5. Future collaborations

    Did you find a bug or a more efficient way for implementing some library's feature?
    Any help for improving the PbDLib are welcome!
    Please contact us or visit http://github.com/ddetommaso/pbdlib

