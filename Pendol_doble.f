	IMPLICIT NONE
	integer N,k,i
	double precision l1,l2,g,m1,m2,t(1000000),pi
	double precision ti,tf,y1i,y2i,v1i,v2i
	double precision y1(1000000),y2(1000000)
	double precision v1(1000000),v2(1000000)
	double precision x1,z1,x2,z2,xa,ya,xb,yb
	external a

	pi=acos(-1.d0)
	N=100000
	l1=1.d0
	l2=l1
	m1=1.d0
	m2=m1
	g=9.81d0
	ti=0.d0
	tf=30.d0
	v1i=0.d0
	v2i=0.d0
	open(1,file="data.dat")

c-cas 1
	write(1,*)"#CAS A: pi/4"	
	y1i=pi/4.d0
	y2i=y1i
	call EULER(N,ti,tf,y1i,v1i,y2i,v2i,
     &			l1,l2,g,a,y1,v1,y2,v2,t)
	do i=1,10000
	 k=10*i
	 xa=x1(y1(k),l1)
	 xb=x2(y1(k),y2(k),l1,l2)
	 ya=z1(y1(k),l1)
	 yb=z2(y1(k),y2(k),l1,l2)
	 write(1,*)xa,ya,xb,yb
	end do
	write(1,*)
	write(1,*)

c-cas 2
	write(1,*)"#CAS B: pi/4+0.1 rads"	
	y1i=0.1d0+pi/4.d0
	y2i=y1i
	call EULER(N,ti,tf,y1i,v1i,y2i,v2i,
     &			l1,l2,g,a,y1,v1,y2,v2,t)
	do i=1,10000
	 k=10*i
	 xa=x1(y1(k),l1)
	 xb=x2(y1(k),y2(k),l1,l2)
	 ya=z1(y1(k),l1)
	 yb=z2(y1(k),y2(k),l1,l2)
	 write(1,*)xa,ya,xb,yb
	end do
	write(1,*)
	write(1,*)

c-cas 2
	write(1,*)"#CAS C: pi/4+0.1 rads"	
	y1i=0.1d0+pi/4.d0
	y2i=y1i+0.1d0
	call EULER(N,ti,tf,y1i,v1i,y2i,v2i,
     &			l1,l2,g,a,y1,v1,y2,v2,t)
	do i=1,10000
	 k=10*(i-1)+1
	 xa=x1(y1(k),l1)
	 xb=x2(y1(k),y2(k),l1,l2)
	 ya=z1(y1(k),l1)
	 yb=z2(y1(k),y2(k),l1,l2)
	 write(1,*)xa,ya,xb,yb
	end do
	write(1,*)
	write(1,*)

	call system("gnuplot pdoble-gif.gnu")
	

	END
c-----------------------------------------------------------------------------------------------------------------------------------------------------------
	subroutine EULER(N,ti,tf,y1i,v1i,y2i,v2i,
     &			l1,l2,g,a,y1,v1,y2,v2,t)
	implicit none
	 integer N,k
	 double precision t(1000000)
	 double precision y1(1000000),v1(1000000)
	 double precision y2(1000000),v2(1000000)
	 double precision h,ti,tf,y1i,v1i,a1,a2,y2i,v2i
	 double precision l1,l2,m1,m2,g
	 h=(tf-ti)/N
	 t(1)=ti
	 y1(1)=y1i
	 v1(1)=v1i
	 y2(1)=y2i
	 v2(1)=v2i
	 call a(l1,l2,m1,m2,g,t(1),y1(1),v1(1),y2(1),v2(1),a1,a2)
	 v1(2)=v1(1)+h*a1
	 y1(2)=y1(1)+h*v1(1)
	 v2(2)=v2(1)+h*a2
	 y2(2)=y2(1)+h*v2(1)
	 do k=2,N+1
	  t(k)=ti+h*(k-1)
	  call a(l1,l2,m1,m2,g,t(k),y1(k),v1(k),y2(k),v2(k),a1,a2)
	  v1(k+1)=v1(k-1)+2.d0*h*a1
	  y1(k+1)=y1(k-1)+2.d0*h*v1(k)
	  v2(k+1)=v2(k-1)+2.d0*h*a2
	  y2(k+1)=y2(k-1)+2.d0*h*v2(k)
	 end do
	return
	end

	subroutine a(l1,l2,m1,m2,g,t,y1,v1,y2,v2,a1,a2)
	implicit none
	 double precision t,y1,v1,y2,v2,a1,a2,l1,l2,g,m1,m2
	 double precision p1,p2,p3,p4,p5,p6,p7,p8,p9
	  p1=-g*(2.d0*m1+m2)*sin(y1)
	  p2=-g*m2*sin(y1-2.d0*y2)
	  p3=-2.d0*sin(y1-y2)*m2
	  p4=l2*(v2**2.d0)+l1*(v1**2.d0)*cos(y1-y2)
	  p5=2.d0*m1+m2-m2*cos(2.d0*(y1-y2))
	 a1=(p1+p2+p3*p4)/(l1*p5)
	  p6=2.d0*sin(y1-y2)
	  p7=l1*(v1**2.d0)*(m1+m2)
	  p8=g*(m1+m2)*cos(y1)
	  p9=l2*m2*(v2**2.d0)*cos(y1-y2)
	 a2=(p6*(p7+p8+p9))/(l2*p5)
	return
	end

	double precision function x1(y1,l1)
	implicit none
	 double precision y1,l1
	 x1=l1*sin(y1)
	return
	end

	double precision function z1(y1,l1)
	implicit none
	 double precision y1,l1
	 z1=-l1*cos(y1)
	return
	end

	double precision function x2(y1,y2,l1,l2)
	implicit none
	 double precision y1,l1,y2,l2
	 x2=l1*sin(y1)+l2*sin(y2)
	return
	end

	double precision function z2(y1,y2,l1,l2)
	implicit none
	 double precision y1,l1,y2,l2
	 z2=-l1*cos(y1)-l2*cos(y2)
	return
	end

	


