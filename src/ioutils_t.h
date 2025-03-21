/*
 *  Copyright Finn Lindgren (2010-2024)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef _FMESH_IOUTILS_T_
#define _FMESH_IOUTILS_T_ 1

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "fmesher_debuglog.h"

namespace fmesh {

// No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R

/*
template <class T>
IOHeader& IOHeader::def(const T& ref)
{
  def();
  return *this;
}
*/

template <class T>
IOHeader &IOHeader::dense(const Matrix<T> &M, IOMatrixtype matrixt) {
  datatype = IODatatype::Dense;
  matrixtype = matrixt;
  switch (matrixt) {
  case IOMatrixtype::Invalid:
    break;
  case IOMatrixtype::General:
    elems = M.rows() * M.cols();
    rows = M.rows();
    cols = M.cols();
    break;
  case IOMatrixtype::Symmetric:
    if (M.rows() <= M.cols()) {
      elems = (M.rows() * (M.rows() + 1)) / 2;
      rows = M.rows();
      cols = M.rows();
    } else {
      elems = (M.cols() * (M.cols() + 1)) / 2;
      rows = M.cols();
      cols = M.cols();
    }
    break;
  case IOMatrixtype::Diagonal:
    if (M.rows() <= M.cols()) {
      elems = M.rows();
      rows = M.rows();
      cols = M.rows();
    } else {
      elems = M.cols();
      rows = M.cols();
      cols = M.cols();
    }
    break;
  }
  return *this;
}

template <class T>
IOHeader &IOHeader::sparse(const SparseMatrix<T> &M, IOMatrixtype matrixt) {
  datatype = IODatatype::Sparse;
  matrixtype = matrixt;
  elems = M.nnz(matrixt);
  switch (matrixt) {
  case IOMatrixtype::Invalid:
    break;
  case IOMatrixtype::General:
    rows = M.rows();
    cols = M.cols();
    break;
  case IOMatrixtype::Symmetric:
    if (M.rows() <= M.cols()) {
      rows = M.cols();
      cols = M.cols();
    } else {
      rows = M.rows();
      cols = M.rows();
    }
    break;
  case IOMatrixtype::Diagonal:
    if (M.rows() <= M.cols()) {
      rows = M.rows();
      cols = M.rows();
    } else {
      rows = M.cols();
      cols = M.cols();
    }
    break;
  }
  return *this;
}

template <class T> IOHelper<T> &IOHelper<T>::OH(std::ostream &output) {
  if (binary_) {
    int header_size = sizeof(IOHeader);
    int header_length = header_size / sizeof(int);
    output.write((char *)&header_length, sizeof(header_length));
    output.write((const char *)&h_, header_size);
  } else {
    output << h_;
    output << std::endl;
  }
  return *this;
}

template <class T> IOHelper<T> &IOHelper<T>::IH(std::istream &input) {
  if (binary_) {
    int header_size = sizeof(IOHeader);
    int file_header_length;
    input.read((char *)&file_header_length, sizeof(file_header_length));
    int file_header_size = file_header_length * sizeof(int);
    if (file_header_size < header_size) {
      h_.dense(Matrix<int>(0), IOMatrixtype::General);
      input.read((char *)&h_, file_header_size);
    } else {
      input.read((char *)&h_, header_size);
      if (file_header_size > header_size) {
        char *buf = new char[file_header_size - header_size];
        input.read(buf, file_header_size - header_size);
        delete[] buf;
      }
    }
  } else {
    input >> h_;
  }
  return *this;
}

template <class T> IOHelper<T> &IOHelper<T>::IH(const IOHeader &h) {
  h_ = h;
  return *this;
}

template <class T> IOHelperM<T> &IOHelperM<T>::OD(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  const bool bin_(IOHelper<T>::is_binary());
  if ((!((h.rows > 0) && (h.cols > 0))) || (!cM_)) {
    return *this;
  }
  if (bin_) {
    switch ((IOMatrixtype)h.matrixtype) {
    case IOMatrixtype::Invalid:
      break;
    case IOMatrixtype::General:
      if ((h.storagetype == IOStoragetype::Rowmajor) || ((*cM_).cols() == 1)) {
        output.write((char *)(*cM_).raw(), sizeof(T) * h.rows * h.cols);
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < h.rows; i++) {
            output.write((char *)&(*cM_)[i][j], sizeof(T));
          }
        }
      }
      break;
    case IOMatrixtype::Symmetric:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          const T *Mrow = (*cM_)[i];
          output.write((char *)&Mrow[i], sizeof(T) * (h.cols - i));
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < j; i++) {
            output.write((char *)&(*cM_)[i][j], sizeof(T));
          }
        }
      }
      break;
    case IOMatrixtype::Diagonal:
      for (int i = 0; i < h.rows; i++) {
        output.write((char *)&(*cM_)[i][i], sizeof(T));
      }
      break;
    }
  } else { /* Text format. */
    output << std::setprecision(15) << std::scientific;
    switch ((IOMatrixtype)h.matrixtype) {
    case IOMatrixtype::Invalid:
      break;
    case IOMatrixtype::General:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          const T *Mrow = (*cM_)[i];
          for (int j = 0; j + 1 < h.cols; j++) {
            output << Mrow[j] << " ";
          }
          output << Mrow[h.cols - 1] << std::endl;
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i + 1 < h.rows; i++) {
            output << (*cM_)[i][j] << " ";
          }
          output << (*cM_)[h.rows - 1][j] << std::endl;
        }
      }
      break;
    case IOMatrixtype::Symmetric:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          const T *Mrow = (*cM_)[i];
          for (int j = i; j + 1 < h.cols; j++) {
            output << Mrow[j] << " ";
          }
          output << Mrow[h.cols - 1] << std::endl;
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < j; i++) {
            output << (*cM_)[i][j] << " ";
          }
          output << (*cM_)[j][j] << std::endl;
        }
      }
      break;
    case IOMatrixtype::Diagonal:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          output << (*cM_)[i][i] << std::endl;
        }
      } else {
        for (int i = 0; i + 1 < h.rows; i++) {
          output << (*cM_)[i][i] << " ";
        }
        output << (*cM_)[h.rows - 1][h.rows - 1] << std::endl;
      }
      break;
    }
  }
  return *this;
}

template <class T> IOHelperM<T> &IOHelperM<T>::ID(std::istream &input) {
  const IOHeader &h(IOHelper<T>::header());
  const bool bin_(IOHelper<T>::is_binary());
  if (!M_) {
    return *this;
  }
  (*M_).clear();
  (*M_).cols(h.cols);
  (*M_).capacity(h.rows);
  if ((h.rows > 0) && (h.cols > 0))
    (*M_)(h.rows - 1, h.cols - 1, T()); /* Initialize last element. */
  if (bin_) {
    switch ((IOMatrixtype)h.matrixtype) {
    case IOMatrixtype::Invalid:
      break;
    case IOMatrixtype::General:
      if ((h.storagetype == IOStoragetype::Rowmajor) || (h.cols == 1)) {
        input.read((char *)(*M_).raw(), sizeof(T) * h.rows * h.cols);
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < h.rows; i++) {
            input.read((char *)&(*M_)(i)[j], sizeof(T));
          }
        }
      }
      break;
    case IOMatrixtype::Symmetric:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          T *Mrow = (*M_)(i);
          input.read((char *)&Mrow[i], sizeof(T) * (h.cols - i));
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < j; i++) {
            input.read((char *)&(*M_)(i)[j], sizeof(T));
          }
        }
      }
      break;
    case IOMatrixtype::Diagonal:
      for (int i = 0; i < h.rows; i++) {
        input.read((char *)&(*M_)(i)[i], sizeof(T));
      }
      break;
    }
  } else { /* Text format. */
    switch ((IOMatrixtype)h.matrixtype) {
    case IOMatrixtype::Invalid:
      break;
    case IOMatrixtype::General:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          T *Mrow = (*M_)(i);
          for (int j = 0; j < h.cols; j++) {
            input >> Mrow[j];
          }
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < h.rows; i++) {
            input >> (*M_)(i)[j];
          }
        }
      }
      break;
    case IOMatrixtype::Symmetric:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          T *Mrow = (*M_)(i);
          for (int j = i; j < h.cols; j++) {
            input >> Mrow[j];
          }
        }
      } else {
        for (int j = 0; j < h.cols; j++) {
          for (int i = 0; i < j + 1; i++) {
            input >> (*M_)(i)[j];
          }
        }
      }
      break;
    case IOMatrixtype::Diagonal:
      if (h.storagetype == IOStoragetype::Rowmajor) {
        for (int i = 0; i < h.rows; i++) {
          input >> (*M_)(i)[i];
        }
      } else {
        for (int i = 0; i < h.rows; i++) {
          input >> (*M_)(i)[i];
        }
      }
      break;
    }
  }
  return *this;
}

template <class T> IOHelperSM<T> &IOHelperSM<T>::OD(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  const bool bin_(IOHelper<T>::is_binary());
  if (!cM_) {
    return *this;
  }
  if (h.storagetype == IOStoragetype::Rowmajor) {
    Matrix1<SparseMatrixTriplet<T>> MT;
    (*cM_).tolist(MT, (IOMatrixtype)h.matrixtype);
    IOHelperM<SparseMatrixTriplet<T>>().cD(&MT).binary(bin_).rowmajor().OD(
        output);
  } else {
    Matrix1int Mr;
    Matrix1int Mc;
    Matrix1<T> Mv;
    (*cM_).tolist(Mr, Mc, Mv, (IOMatrixtype)h.matrixtype);
    IOHelperM<int>().cD(&Mr).binary(bin_).colmajor().OD(output);
    IOHelperM<int>().cD(&Mc).binary(bin_).colmajor().OD(output);
    IOHelperM<T>().cD(&Mv).binary(bin_).colmajor().OD(output);
  }
  return *this;
}

template <class T> IOHelperSM<T> &IOHelperSM<T>::ID(std::istream &input) {
  const IOHeader &h(IOHelper<T>::header());
  const bool bin_(IOHelper<T>::is_binary());
  if (!M_) {
    return *this;
  }
  (*M_).clear();
  if (h.elems == 0) {
    return *this;
  }
  if (h.storagetype == IOStoragetype::Rowmajor) {
    Matrix1<SparseMatrixTriplet<T>> MT;
    MT(h.elems - 1) = SparseMatrixTriplet<T>();
    IOHelperM<SparseMatrixTriplet<T>>().D(&MT).binary(bin_).rowmajor().ID(
        input);
    (*M_).fromlist(MT, (IOMatrixtype)h.matrixtype);
  } else {
    Matrix1int Mr;
    Mr(h.elems - 1) = 0;
    Matrix1int Mc;
    Mc(h.elems - 1) = 0;
    Matrix1<T> Mv;
    Mv(h.elems - 1) = T();
    IOHelperM<int>().D(&Mr).binary(bin_).colmajor().ID(input);
    IOHelperM<int>().D(&Mc).binary(bin_).colmajor().ID(input);
    IOHelperM<T>().D(&Mv).binary(bin_).colmajor().ID(input);
    (*M_).fromlist(Mr, Mc, Mv, (IOMatrixtype)h.matrixtype);
  }
  return *this;
}

template <class T> IOHelperM<T> &IOHelperM<T>::OH_2009(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  output << h.rows;
  output << std::endl;
  return *this;
}

template <class T> IOHelperSM<T> &IOHelperSM<T>::OH_2009(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  if ((IOMatrixtype)h.matrixtype == IOMatrixtype::Diagonal) {
    output << h.rows;
    output << std::endl;
  }
  return *this;
}

template <class T> IOHelperM<T> &IOHelperM<T>::OD_2009(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  if ((!((h.rows > 0) && (h.cols > 0))) || (!cM_)) {
    return *this;
  }
  output << std::setprecision(15) << std::scientific;
  for (int i = 0; i < h.rows; i++) {
    const T *Mrow = (*cM_)[i];
    output << i << " ";
    for (int j = 0; j + 1 < h.cols; j++) {
      output << Mrow[j] << " ";
    }
    output << Mrow[h.cols - 1] << std::endl;
  }
  return *this;
}

template <class T> IOHelperSM<T> &IOHelperSM<T>::OD_2009(std::ostream &output) {
  const IOHeader &h(IOHelper<T>::header());
  if (!cM_) {
    return *this;
  }
  if ((IOMatrixtype)h.matrixtype == IOMatrixtype::Diagonal) {
    Matrix1<T> MT;
    for (int r = 0; r < (*cM_).rows(); r++)
      MT(r) = (*cM_)[r][r];
    IOHelperM<T>().cD(&MT).OD_2009(output);
  } else {
    Matrix1<SparseMatrixTriplet<T>> MT;
    (*cM_).tolist(MT, (IOMatrixtype)h.matrixtype);
    IOHelperM<SparseMatrixTriplet<T>>().cD(&MT).ascii().rowmajor().OD(output);
  }
  return *this;
}

template <class T> IOHeader &IOHeader::def(const T &ref) {
  (void)(ref);
  def();
  valuetype = IOValuetype::Invalid;
  return *this;
}

template <class T> IOHeader::IOHeader(const T &ref) { def(ref); }

template <class T>
void save_M(std::string filename, const Matrix<T> &M, MCCInfo mccinfo,
            bool binary) {
  std::ofstream O;
  O.open(filename.c_str(),
         (binary ? (std::ios::out | std::ios::binary) : std::ios::out));
  IOHelperM<T> ioh;
  ioh.cD(&M).matrixtype((IOMatrixtype)mccinfo.matrixtype);
  ioh.binary(binary).OH(O).OD(O);
  O.close();
}

template <class T>
void save_SM(std::string filename, const SparseMatrix<T> &M, MCCInfo mccinfo,
             bool binary) {
  std::ofstream O;
  O.open(filename.c_str(),
         (binary ? (std::ios::out | std::ios::binary) : std::ios::out));
  IOHelperSM<T> ioh;
  ioh.cD(&M).matrixtype((IOMatrixtype)mccinfo.matrixtype);
  ioh.binary(binary).OH(O).OD(O);
  O.close();
}

template <class T>
void MatrixC::input_raw_M(std::istream &input, Matrix<T> &M) const {
  M.load_ascii_2009(input);
}

#endif // not FMESHER_WITH_R


} /* namespace fmesh */

#ifdef FMESHER_WITH_R
namespace Rcpp {

template<>
inline SEXP wrap(const fmesh::MatrixC& obj) {
  return obj.Rcpp_wrap();
}

}
#endif

#endif
