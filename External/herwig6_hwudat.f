      subroutine herwig_init
      integer lnhwrt,lnhrd,lnhout,lnhdcy
      common/heplun/lnhwrt,lnhrd,lnhout,lnhdcy
      external hwudat
      integer n
      integer istr,nevt

      lnhwrt=0
      lnhrd=0
      lnhdcy=0
      lnhout=22
      lhwout=lnhout
      return
      end
