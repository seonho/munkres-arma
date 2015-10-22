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

#if !defined(_MUNKRES_HPP_)
#define _MUNKRES_HPP_

#pragma once

#include <armadillo>
#include <cassert>
#include <list>

namespace arma
{
	//! Modulus after division
	template <typename vec_type>
	inline vec_type mod(const vec_type& X, typename vec_type::elem_type Y)
	{
		typedef typename vec_type::elem_type elem_type;
		assert(Y != 0);
		
		vec_type M;
		switch (X.vec_state) {
		case 0: // matrix
			M = X - arma::conv_to<vec_type>::from(
				arma::floor(arma::conv_to<mat>::from(X) / (double)Y)) * Y;
			break;
		case 1:
		case 2:
			M = X - arma::conv_to<vec_type>::from(
				arma::floor(arma::conv_to<vec>::from(X) / (double)Y)) * Y;
			break;
		}
		
		return M;
	}

	arma::umat ind2sub(const arma::SizeMat& siz, const arma::uvec& ndx)
	{
		arma::umat sub(ndx.size(), 2);
		sub.col(0) = mod(ndx, siz.n_rows);
		sub.col(1) = ndx / siz.n_rows;
		return sub;
	}

	template <typename mat_type>
	inline mat_type sortrows(const mat_type& X, uword col)
	{
		return X.rows(arma::stable_sort_index(X.col(col)));
	}
};

template <typename eT>
class munkres
{
	typedef arma::uword		size_type;
	typedef arma::uchar_mat	ucmat;

	static const unsigned char NORMAL = 0;
	static const unsigned char STAR	= 1;
	static const unsigned char PRIME	= 2;

public:
	arma::umat solve(const arma::Mat<eT>& m)
	{
		const size_type rows = m.n_rows, 
			columns = m.n_cols, 
			size = std::max(rows, columns);

		matrix_ = m;

		if (!matrix_.is_square()) {
			const eT value = matrix_.max();
			matrix_.resize(size, size);

			if (rows > columns)
				matrix_.cols(columns, size - 1).fill(value);
			else
				matrix_.rows(rows, size - 1).fill(value);
		}

		mask_.set_size(size, size);
		mask_.fill(NORMAL);

		row_mask_.clear();
		row_mask_.resize(size);

		col_mask_.clear();
		col_mask_.resize(size);

		// Prepare the matrix values...

        // If there were any infinities, replace them with a value greater
        // than the maximum value in the matrix.
		if (matrix_.has_inf())
			matrix_.elem(find_nonfinite(matrix_)).fill(matrix_.elem(arma::find_finite(matrix_)).max());

		// minimize along direction
		matrix_.each_col() -= min(matrix_, 1);
		matrix_.each_row() -= min(matrix_, 0);

		// Follow the steps
        int step = 1;
        while ( step ) {
            switch ( step ) {
            case 1:
                step = step1();
                // step is always 2
                break;
            case 2:
                step = step2();
                // step is always either 0 or 3
                break;
            case 3:
                step = step3();
                // step in [3, 4, 5]
                break;
            case 4:
                step = step4();
                // step is always 2
                break;
            case 5:
                step = step5();
                // step is always 3
                break;
            }
        }
		
		// Remove the excess rows or columns that we added to fit the
        // input to a square matrix.
		matrix_.resize(rows, columns);
		mask_.resize(rows, columns);

		// sort subscripts by row subscript.
		return arma::sortrows(arma::ind2sub(arma::size(mask_), find(mask_ == STAR)), 0);
	}

private:
	int step1()
	{
		const size_type rows = matrix_.n_rows,
			columns = matrix_.n_cols;

		for ( size_type row = 0 ; row < rows ; row++ ) {
			for ( size_type col = 0 ; col < columns ; col++ ) {
				if ( 0 == matrix_.at(row, col) ) {
					bool starred = false;
					for ( size_type nrow = 0 ; nrow < rows ; nrow++ ) {
						if ( STAR == mask_.at(nrow, col) ) {
							starred = true;
							break;
						}
					}

					if (!starred) {
						for ( size_type ncol = 0 ; ncol < columns ; ncol++ ) {
							if ( STAR == mask_.at(row, ncol) ) {
								starred = true;
								break;
							}
						}
					}

					if (!starred)
						mask_.at(row, col) = STAR;
				}
			}
		}

		return 2;
	}

	int step2()
	{
		const size_type rows = matrix_.n_rows,
			columns = matrix_.n_cols;
		size_type covercount = 0;

		for ( size_type row = 0 ; row < rows ; row++ ) {
			for ( size_type col = 0 ; col < columns ; col++ ) {
				if ( STAR == mask_.at(row, col) ) {
					col_mask_[col] = true;
					covercount++;
				}
			}
		}

		
		if ( covercount >= min(size(matrix_)) ) {
#ifdef _DEBUG
			std::cout << "Final cover count: " << covercount << std::endl;
#endif
			return 0;
		}

#ifdef _DEBUG
		std::cout << "Munkres matrix has " << covercount << " of " << min(size(matrix_)) << " columns covered:" << std::endl;
		std::cout << matrix_ << std::endl;
#endif

		return 3;
	}

	int step3()
	{
		/*
		Main Zero Search
		1. Find an uncovered Z in the distance matrix and prime it. If no such zero exists, go to Step 5
		2. If No Z* exists in the row of the Z', go to Step 4.
		3. If a Z* exists, cover this row and uncover the column of the Z*. Return to Step 3.1 to find a new Z
		*/

		if (find_uncovered(0, saverow_, savecol_))
			mask_.at(saverow_, savecol_) = PRIME; // prime it.
		else
			return 5;

		for ( size_type col = 0 ; col < matrix_.n_cols ; col++ ) {
			if ( mask_.at(saverow_, col) == STAR ) {
				row_mask_[saverow_]	= true;  //cover this row and
				col_mask_[col]		= false; // uncover the column containing the starred zero

				return 3; // repeat
			}
		}

		return 4; // no starred zero in the row containing this primed zero
	}

	int step4()
	{
		const size_type rows = matrix_.n_rows,
			columns = matrix_.n_cols;

		// seq contains pairs of row/column values where we have found
		// either a star or a prime that is part of the ``alternating sequence``.
		std::list<std::pair<int, int> > seq;
		// use saverow, savecol from step 3.
		std::pair<int, int> z0(saverow_, savecol_);
		seq.insert(seq.end(), z0);

		// We have to find these two pairs:
		std::pair<int, int> z1(-1, -1);
		std::pair<int, int> z2n(-1, -1);

		size_type row, col = savecol_;
		/*
		Increment Set of Starred Zeros

		1. Construct the ``alternating sequence'' of primed and starred zeros:

			Z0 : Unpaired Z' from Step 4.2
			Z1 : The Z* in the column of Z0
			Z[2N] : The Z' in the row of Z[2N-1], if such a zero exists
			Z[2N+1] : The Z* in the column of Z[2N]

		The sequence eventually terminates with an unpaired Z' = Z[2N] for some N.
		*/
		bool madepair;
		do {
			madepair = false;
			for ( row = 0 ; row < rows ; row++ ) {
				if ( mask_.at(row, col) == STAR ) {
					z1.first = row;
					z1.second = col;
					
					if ( pair_in_list(z1, seq) )
						continue;

					madepair = true;
					seq.insert(seq.end(), z1);
					break;
				}
			}

			if ( !madepair )
				break;

			madepair = false;

			for ( col = 0 ; col < columns ; col++ ) {
				if ( mask_.at(row, col) == PRIME ) {
					z2n.first = row;
					z2n.second = col;

					if ( pair_in_list(z2n, seq) )
						continue;

					madepair = true;
					seq.insert(seq.end(), z2n);
					break;
				}
			}
		} while ( madepair );

		std::for_each(seq.begin(), seq.end(), [&](std::pair<int, int>& p) {
			// 2. Unstar each starred zero of the sequence.
			if ( mask_.at(p.first, p.second) == STAR )
				mask_.at(p.first, p.second) = NORMAL;

			// 3. Star each primed zero of the sequence,
			// thus increasing the number of starred zeros by one.
			if ( mask_.at(p.first, p.second) == PRIME )
				mask_.at(p.first, p.second) = STAR;
		});

		// 4. Erase all primes, uncover all columns and rows,
		mask_.elem(arma::find(mask_ == PRIME)).fill(NORMAL);

		std::fill(row_mask_.begin(), row_mask_.end(), false);
		std::fill(col_mask_.begin(), col_mask_.end(), false);

		// and return to Step 2.
		return 2;
	}

	int step5()
	{
		const size_type rows = matrix_.n_rows,
			columns = matrix_.n_cols;

		/*
		New Zero Manufactures

		1. Let h be the smallest uncovered entry in the (modified) distance matrix_.
		2. Add h to all covered rows.
		3. Subtract h from all uncovered columns
		4. Return to Step 3, without altering stars, primes, or covers.
		*/

		eT h = 0;

		for ( size_type row = 0 ; row < rows ; row++ ) {
			if ( !row_mask_[row] ) {
				for ( size_type col = 0 ; col < columns ; col++ ) {
					if ( !col_mask_[col] ) {
						if ( (h > matrix_.at(row, col) && matrix_.at(row, col) != 0) || h == 0 )
							h = matrix_.at(row, col);
					}
				}
			}
		}

		for ( size_type row = 0 ; row < rows ; row++ ) {
			if ( row_mask_[row] )
				matrix_.row(row) += h;
		}

		for ( size_type col = 0 ; col < columns ; col++ ) {
			if ( !col_mask_[col] )
				matrix_.col(col) -= h;
		}

		return 3;
	}

	inline bool pair_in_list(const std::pair<int,int> &needle, const std::list<std::pair<int,int> > &haystack) const
	{
		return std::find(haystack.begin(), haystack.end(), needle) != haystack.end();
	}

	inline bool find_uncovered(const eT item, size_type& row, size_type& col) const
	{
		const size_type rows = matrix_.n_rows,
			columns = matrix_.n_cols;

		// adjusted to column major ordering
		for ( col = 0 ; col < columns ; col++ ) {
			if ( !col_mask_[col] ) {
				for ( row = 0 ; row < rows ; row++ ) {
					if ( !row_mask_[row] ) {
						if ( matrix_.at(row, col) == item )
							return true;
					}
				}
			}
		}

		return false;
	}

	arma::Mat<eT>		matrix_;
	ucmat				mask_;
	
	std::vector<bool>	row_mask_;
	std::vector<bool>	col_mask_;

	size_type			saverow_;
	size_type			savecol_;
};

#endif /* !defined(_MUNKRES_HPP_) */