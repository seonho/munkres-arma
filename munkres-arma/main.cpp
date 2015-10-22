/*
 *   Copyright (c) 2007 John Weaver
 *   Copyright (c) 2015 Seonho Oh
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#include "munkres.hpp"

int main(int argc, char* argv[])
{
	int nrows = 4,
		ncols = 3;

	if ( argc == 3 ) {
		nrows = atoi(argv[1]);
		ncols = atoi(argv[2]);
	}

	// Initialize matrix with random values.
	arma::mat cost = arma::conv_to<arma::mat>::from(arma::randi(nrows, ncols, arma::distr_param(1, 50)));
	cost.print("cost = ");

	// Apply Munkres algorithm to cost matrix.
	munkres<double> m;
	auto assignments = m.solve(cost);
	assignments.print("assignments = ");

	return 0;
}