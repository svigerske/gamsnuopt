      subroutine attrvr(i,idir,dpc,upc,prio,iunt)
      integer i,idir
     #     ,iunt
      double precision dpc,upc,prio
c
      call attrvrbc(i,idir,dpc,upc,prio,iunt)
c     
      return
      end

      subroutine deffun(pbnam,nfmax,nfnc,kfun,cl,cu,funam,khess,ir)
      integer nfmax,nfnc,kfun(nfmax),khess(nfmax),ir
      double precision cl(nfmax),cu(nfmax)
      character*40 pbnam
      character*8 funam(nfmax)
c
      integer i
c
      call deffunbc(pbnam,nfmax,nfnc,kfun,cl,cu,funam,khess,ir)
c
      return
      end

      subroutine defsdp(nblkmx,nblk,dim,clsdp,sdpnam
     #     ,msdpmx,msdp
     #     ,isdp,jfun,imrow,jmcol,ir)
      integer nblkmx,nblk,dim(nblkmx)
     #     ,msdpmx,msdp
     #     ,isdp(msdpmx),jfun(msdpmx)
     #     ,imrow(msdpmx),jmcol(msdpmx)
     #     ,ir
      double precision clsdp(nblkmx)
      character*8 sdpnam(nblkmx)
c
      nblk = 0
c
      return
      end

      subroutine defvar(nvmax,nvar,kvar,bl,bu,vrnam,ir)
      integer nvmax,nvar,kvar(nvmax),ir
      double precision bl(nvmax),bu(nvmax)
      character*8 vrnam(nvmax)
c     
      integer i
c
      call defvarbc(nvmax,nvar,kvar,bl,bu,vrnam,ir)
c
      return
      end

      subroutine funl(nemax,ne,ifun,jvar,a,ir)
      integer nemax,ne,ifun(nemax),jvar(nemax),ir
      double precision a(nemax)
c
      call funlbc(nemax,ne,ifun,jvar,a,ir)
c
      return
      end

      subroutine funnl(nvar,nfnc,x,f,ir)
      integer nvar,nfnc,ir
      double precision x(nvar),f(nfnc)
c
      call funnlbc(nvar,nfnc,x,f,ir)
c
      return
      end

      subroutine funq(nemax,ne,ifun,j1var,j2var,h,ir)
      integer nemax,ne,ifun(nemax),j1var(nemax),j2var(nemax),ir
      double precision h(nemax)
c
      call funqbc(nemax,ne,ifun,j1var,j2var,h,ir)
c
      return
      end

      subroutine gradnl(nvar,x,nemax,ne,ifun,jvar,a,ir)
      integer nvar,nemax,ne,ifun(nemax),jvar(nemax),ir
      double precision x(nvar),a(nemax)
c
      call gradnlbc(nvar,x,nemax,ne,ifun,jvar,a,ir)
c
      return
      end

      subroutine hessnl(nvar,nfnc,x,y,nemax,ne,j1var,j2var,h,ir)
      integer nvar,nfnc,nemax,ne,j1var(nemax),j2var(nemax),ir
      double precision x(nvar),y(nfnc),h(nemax)
c
      call hessnlbc(nvar,nfnc,x,y,nemax,ne,j1var,j2var,h,ir)
c
      return
      end

      subroutine initvr(nvar,xini,ir)
      integer nvar,ir
      double precision xini(nvar)
c
      call initvrbc(nvar,xini,ir)
c
      return
      end

      subroutine typevr(nvar,ivtyp)
      integer nvar,ivtyp(nvar)
c
c     変数のタイプを与えるインタフェース
c
      call typevrbc(nvar,ivtyp)
c
      return
      end
