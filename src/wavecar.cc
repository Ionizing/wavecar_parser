#include <wavecar.hpp>

// #define DEBUG

namespace ionizing {
  WAVECAR::WAVECAR(const char* FileName) : io(FileName) {
    _info._fileSize = io.getFileSize();

#ifdef DEBUG
    std::cout << "Open WAVECAR: " << FileName << " successfully!\n"
      << "File size == " << _fileSize << std::endl;
#endif

    /*
     * std::cout << "WAVECAR created successfully!" << std::endl;
     */

    io.seek(0);

    read_info();
    read_header();
    read_band();
  }

  void WAVECAR::read_info() {
    _iRec = 0;
    io.seek(0);

    _info._recordLength = static_cast<int>( io.readElement<double>());
    _info._nSpin        = static_cast<int>( io.readElement<double>());
    _info._precisionTag = static_cast<int>( io.readElement<double>());

    if (_info._recordLength <= 0) {
      std::cout << "****ERROR**** WAVECAR _recordLength wrong: "
        << _info._recordLength << " !\n";
      std::abort();
    }
    
    if (_info._nSpin != 0 and _info._nSpin != 1) {
      std::cout << "****ERROR**** WAVECAR _nSpin wrong: "
        << _info._nSpin << " !\n";
      std::abort();
    }

    if (_info._precisionTag <= 0) {
      std::cout << "****ERROR**** WAVECAR _precisionTag wrong: "
        << _info._precisionTag << " !\n";
      std::abort();
    }

    if (_info._isDoubleType) {
      _maxOfNPlaneWaves = (_info._recordLength * 2 - 1) / sizeof(std::complex<double>) / 2;
    } else {
      _maxOfNPlaneWaves = (_info._recordLength * 2 - 1) / sizeof(std::complex< float>) / 2;
    }

#ifdef DEBUG
    std::cout << "\t_maxOfNPlaneWaves == " << _maxOfNPlaneWaves << " .\n";
#endif


    _info._isDoubleType = getPrecisionType();


#ifdef DEBUG
    std::cout << "\t_recordLength     == " << _info._recordLength << " \n"
              << "\t_nSpin            == " << _info._nSpin << " \n"
              << "\t_precisionTag     == " << _info._precisionTag << " \n";
#endif

  } // end of function read_info

  bool WAVECAR::getPrecisionType() const {
    switch (_info._precisionTag) {
      case 45200: return false;
      case 45210: return true;
      case 53310: ;
      case 53300: std::cerr << "This is a VASP5 WAVECAR, not implemented yet\n";
                  std::abort();
      default   : std::cerr << "Invalid WAVECAR precision tag: " << _info._precisionTag << " !\n";
                  std::abort();
    }
  }

  void WAVECAR::read_header() {
    _iRec = 1;
    io.seek(_iRec * _info._recordLength);
    _header._nKpoints = static_cast<int>( io.readElement<double>());
    _header._nBands   = static_cast<int>( io.readElement<double>());
    _header._enCut    = io.readElement<double>();

    _header._latticeVectors    = io.readMatrix<double>(3, 3);
    _header._reciprocalVectors = _header._latticeVectors.inverse();
    _header._omega             = _header._latticeVectors.determinant();
    _header._eFermi            = io.readElement<double>();

    // We cannot check everything for the code style and performance issues.
    // so I only check the validity of nKpoints and nBands;
    if (_header._nKpoints <= 0 or _header._nBands <= 0) {
      std::cout << "****ERROR**** _nKpoints or nBands error: "
        << _header._nKpoints << "\t" << _header._nBands << " !\n";
      std::abort();
    } 

#ifdef DEBUG
    std::cout << "\t_nKpoints         == "   << _header._nKpoints << " \n"
              << "\t_nBands           == "   << _header._nBands << " \n"
              << "\t_enCut            == "   << _header._enCut << " \n"
              << "\t_latticeVectors   == \n" << _header._latticeVectors << " \n"
              << "\t_eFermi           == "   << _header._eFermi << " \n";
#endif
  }

  void WAVECAR::read_band() {
    _iRec = 2;
    _nPlaneWaves  .resize(_header._nKpoints * _info._nSpin);
    _bands        .resize(_info._nSpin);
    _fermiWeights .resize(_info._nSpin);
    _complexWaves .resize(_info._nSpin);
    _kVectors     .resize(_header._nKpoints, 3);
        // std::abort();

    for (int i=0; i!=_info._nSpin; ++i) {
      _bands(i)       .resize(_header._nKpoints,  _header._nBands);         
      _fermiWeights(i).resize(_header._nKpoints,  _header._nBands);   
      _complexWaves(i).resize(_header._nKpoints, _maxOfNPlaneWaves);  
    } // Initialize _bands, _fermiWeights and _complexWaves


    for (int ispin=0; ispin!=_info._nSpin; ++ispin) {
      for (int kpoint=0; kpoint!=_header._nKpoints; ++kpoint) {
        io.seek(_iRec * _info._recordLength);
        VectorXd buf1 = io.readVectorRow<double>(4);  


#ifdef DEBUG
        std::cout << "\t\tbuf1 == " << buf1.transpose() << std::endl;
#endif


        Matrix<double, Dynamic, 3> buf2 = 
            io.readMatrix<double>(_header._nBands, 3);  


#ifdef DEBUG
        std::cout << "\n\t buf2 = \n" << buf2 << std::endl;
#endif


        // buf2 is stored as:   Re(band[0]), Im(band[0]), fermiWeights[0];
        //  in RowMajor style   Re(band[1]), Im(band[1]), fermiWeights[1];
        //                      ...           ...         ...

        // std::abort();

        int current_nplw = 0;
        if (0 == ispin) {
          current_nplw          = static_cast<int>(buf1(0));
          _nPlaneWaves(kpoint)  = current_nplw;
          _kVectors.row(kpoint) = buf1.segment<3>(1);


#ifdef DEBUG
          std::cout << _kVectors.row(kpoint) << std::endl;
          std::cout << "\tcurrent_nplw == " << current_nplw << " \n";
#endif


        }

        _bands(ispin).row(kpoint).real() = buf2.col(0); 
        _bands(ispin).row(kpoint).imag() = buf2.col(1); 
        _fermiWeights(ispin).row(kpoint) = buf2.col(2); 

        // Next record is complex CW(:)
        ++_iRec;  // head to complex wavefunction
        for (int iband=0; iband!=_header._nBands; ++iband) {
          io.seek(_iRec * _info._recordLength);
          if (true == _info._isDoubleType) {
            _complexWaves(ispin)(kpoint, iband) = 
                io.readVectorCol<std::complex<double>>(current_nplw);
          } else {
            _complexWaves(ispin)(kpoint, iband) = 
                io.readVectorCol<std::complex< float>>(current_nplw)
                  .cast<std::complex<double>>();
          }


#ifdef DEBUG
          std::cout << "\n\t\t wavefunction for spin " << ispin 
            << " , kpoint " << kpoint <<" , band " << iband << " is:\n" << std::endl;
          std::cout << _complexWaves(ispin)(kpoint, iband) << std::endl;
          if (true == _isDoubleType) {
            std::cout << "DEBUG: test if reached end of planewaves,\n" 
              << io.readVectorCol<std::complex<double>>(3)
              << "\n\n\n";
          } else {
            std::cout << "DEBUG: test if reached end of planewaves,\n" 
              << io.readVectorCol<std::complex< float>>(3)
              << "\n\n\n";
          }
#endif


          ++_iRec;
        } // end iband for of _complexWaves

      } // kpoint for
    } // ispin for

  } // end of read_band();


  // returns recl, nspin, prectag, isdoubletype, filesize
  const WAVECAR::Info& WAVECAR::getInfo() const {
    return _info;
  }

  // returns header: nkpts, nbands, encut, Acell
  const WAVECAR::Header& WAVECAR::getHeader() const {
    return _header;
  }

  // returns lattice vectors in real space
  const Mat33d& WAVECAR::getLatticeVectors() const {
    return _header._latticeVectors;
  } 

  // returns kVectors
  const MatX3d& WAVECAR::getKVectors() const {
    return _kVectors;
  }

  // returns bands datastructure 
  const ColVecT<Matcd>& WAVECAR::getBands() const {
    return _bands;
  }

  // returns fermiweight 3d array
  const ColVecT< Matd>& WAVECAR::getFermiWeights() const {
    return _fermiWeights;
  }

  // returns planewave coeffs
  const ColVecT<MatT<Veccd>>& WAVECAR::getComplexWaves() const {
    return _complexWaves;
  }

  const Mat33d& WAVECAR::getReciprocalVectors() const {
    return _header._reciprocalVectors;
  }

  const VectorXi& WAVECAR::getNPlaneWaves() const {
    return _nPlaneWaves;
  }

  const bool WAVECAR::checkIndex( const int ispin,
         const int ikpoint, const int iband) const {
    if (ispin < 0 or ispin >= 2) {
      std::cerr << "Invalid spin index: " << ispin << std::endl;
      return false;
    }

    if (ikpoint < 0 or ikpoint >= _header._nKpoints) {
      std::cerr << "Invalid kpoint index: " << ikpoint << std::endl;
      return false;
    }

    if (iband < 0 or iband >= _header._nBands) {
      std::cerr << "Invalid band index: " << iband << std::endl;
      return false;
    }

    return true;
  }

  const Veccd WAVECAR::getBandCoeff( const int ispin, 
                                     const int ikpoint,
                                     const int iband, 
                                     const bool is_normed) const {
    if (_complexWaves.size() == 0) {
      std::cerr << __FUNCTION__ << ": _complexWaves is empty!" << std::endl;
      std::abort();
    }

    auto out = _complexWaves(ispin)(ikpoint, iband);

    if (is_normed) {
      /*
       * std::cout << __FILE__ << __LINE__ << "out.norm = " << out.norm() << std::endl;
       */
      out /= out.norm();
    }

    return out;
  }

  void WAVECAR::printInfo(std::ostream& os) const {
    using namespace std;
    // backup initial format flag
    ios_base::fmtflags fmtflag_backup(os.flags());

    os << endl
       << "*****************************************************************************\n"
       << "*********************** WAVECAR INFO FOR CURRENT FILE ***********************\n"
       << "*****************************************************************************\n"
       << setiosflags(ios::scientific | ios::right | ios::showpos | ios::showpoint) 
       << setprecision(10) << boolalpha;

// Info
    os << setw(20) << " Record Length = " << setw(20) << _info._recordLength << endl
       << setw(20) << "         NSpin = " << setw(20) << _info._nSpin        << endl
       << setw(20) << " Precision Tag = " << setw(20) << _info._precisionTag << endl
       << setw(20) << "Is Double Type = " << setw(20) << _info._isDoubleType << endl
       << setw(20) << "     File Size = " << setw(20) << _info._fileSize     << endl ;
// Header 
    os << setw(20) << "     N KPoints = " << setw(20) << _header._nKpoints          << endl
       << setw(20) << "       N Bands = " << setw(20) << _header._nBands            << endl
       << setw(20) << " Energy Cutoff = " << setw(20) << _header._enCut             << endl
       << setw(20) << "       E-Fermi = " << setw(20) << _header._eFermi            << endl
       << setw(20) << "Lattice Vector = " << endl << _header._latticeVectors.format(HeavyFmt)    << endl
       << setw(20) << "Lattice Volume = " << endl << _header._reciprocalVectors.format(HeavyFmt) << endl ;
    // restore the format flag
    os.flags(fmtflag_backup);
    os << "*****************************************************************************\n"
       << "**************************** END OF WAVECAR INFO ****************************\n"
       << "*****************************************************************************\n"
       << endl;
  }

}; // namespace ionizing
