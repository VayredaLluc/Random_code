program double_pendulum
    implicit none
    real(8), parameter :: g=10.d0,l1=1.d0,l2=1.d0,m1=1.d0,m2=1.d0
    real(8), parameter :: dt=0.01d0
    real(8), parameter :: pi=4*atan(1.d0)
    real(8) :: a1,a2,w1,w2
    real(8) :: tt
    real(8) :: xx1,xx2,yy1,yy2,vv1,vv2,EE
    real(8), dimension(4) :: yy
    integer :: i,j
    character(len=30) :: filename

    do i=1,10
        write(filename, '("double_pendulum_res_", I0, ".dat")') i
        open(i,file=filename,status='replace')
        a1=pi/6.d0
        a2=3*pi/4.d0+(i-1)/25.d0
        w1=0.d0
        w2=0.d0
        yy = [a1,a2,w1,w2]
        tt = 0.d0
        do j=1,10000
            xx1 = x1(yy(1),yy(2))
            xx2 = x2(yy(1),yy(2))
            yy1 = y1(yy(1),yy(2))
            yy2 = y2(yy(1),yy(2))
            vv1 = v1(yy(3),yy(4),yy(1),yy(2))
            vv2 = v2(yy(3),yy(4),yy(1),yy(2))
            EE=E(vv1,vv2,yy1,yy2)
            write(i,*) tt,xx1,yy1,xx2,yy2,EE

            call RK4_step(der_dble_pend,tt,yy,dt)
        end do
        close(i)
    end do






contains
    ! time derivatives (first angular velocities and then angular accelerations)
    function d_a1(w1) result(w)
        implicit none
        real(8) :: w1,w
        w = w1
    end function d_a1

    function d_a2(w2) result(w)
        implicit none
        real(8) :: w2,w
        w = w2
    end function d_a2

    function d_w1(a1,a2,w1,w2) result(acc)
        implicit none
        real(8) :: a1,a2,w1,w2
        real(8) :: delta,num,den,acc
        delta = a1-a2
        num = m2*l1*w1**2.d0*sin(2.d0*delta)/2.d0+m2*l2*w2**2.d0*sin(delta)+m1*g*sin(a1)+m2*g*sin(delta)*cos(a2)
        den = -l1*(m1+m2*sin(delta)**2.d0)
        acc = num/den
    end function d_w1

    function d_w2(a1,a2,w1,w2) result(acc)
        implicit none
        real(8) :: a1,a2,w1,w2
        real(8) :: delta,num,den,acc
        delta = a1-a2
        num = m2*l2*w2**2.d0*sin(2.d0*delta)/2.d0+(m1+m2)*l1*w1**2.d0*sin(delta)+(m1+m2)*g*sin(delta)*cos(a1)
        den = l2*(m1+m2*sin(delta)**2.d0)
        acc = num/den
    end function d_w2

    ! Function with all the derivatives of the double pendulum (the equations of motion)
    function der_dble_pend(t,y) result(d_y)
        implicit none
        real(8) :: a1,a2,w1,w2
        real(8) :: t
        real(8), dimension(:) :: y
        real(8), dimension(size(y)) :: d_y
        a1 = y(1)
        a2 = y(2)
        w1 = y(3)
        w2 = y(4)

        d_y(1)=d_a1(w1)
        d_y(2)=d_a2(w2)
        d_y(3)=d_w1(a1,a2,w1,w2)
        d_y(4)=d_w2(a1,a2,w1,w2)
    end function der_dble_pend
        

    ! Runge-Kutta method
    subroutine RK4_step(der_y,t,y,h)
        implicit none

        interface
            function der_y(time, vector) result(der_vector)
                implicit none
                real(8) :: time
                real(8), dimension(:) :: vector
                real(8), dimension(size(vector)) :: der_vector
            end function der_y
        end interface
        
        real(8), dimension(:) :: y
        real(8) :: t,h
        real(8), dimension(size(y)) :: k1,k2,k3,k4

        k1 = der_y(t,y)
        k2 = der_y(t+0.5d0*h,y+0.5d0*h*k1)
        k3 = der_y(t+0.5d0*h,y+0.5d0*h*k2)
        k4 = der_y(t+h,y+h*k3)

        t = t+h
        y = y + h*(k1+2.d0*k2+2.d0*k3+k4)/6.d0

    end subroutine

    ! cartesian coordinates
    function x1(a1,a2) result(x)
        implicit none
        real(8) :: a1,a2,x
        x = l1*sin(a1)
    end function x1

    function y1(a1,a2) result(y)
        implicit none
        real(8) :: a1,a2,y
        y = -l1*cos(a1)
    end function y1

    function x2(a1,a2) result(x)
        implicit none
        real(8) :: a1,a2,x
        x = l1*sin(a1)+l2*sin(a2)
    end function x2

    function y2(a1,a2) result(y)
        implicit none
        real(8) :: a1,a2,y
        y = -l1*cos(a1)-l2*cos(a2)
    end function y2

    !velocities, momentum and energy
    function v1(w1,w2,a1,a2) result(vel)
        implicit none
        real(8) :: w1,w2,a1,a2,vel
        vel = ((l1*w1)**2.d0)**0.5d0
    end function v1

    function v2(w1,w2,a1,a2) result(vel)
        implicit none
        real(8) :: w1,w2,a1,a2,vel
        vel = ((l1*w1)**2.d0+(l2*w2)**2.d0+2.d0*l1*l2*w1*w2*cos(a1-a2))**0.5d0
    end function v2

    function p1(vel) result(mom)
        implicit none
        real(8) :: vel,mom
        mom = vel*m1
    end function p1

    function p2(vel) result(mom)
        implicit none
        real(8) :: vel,mom
        mom = vel*m2
    end function p2

    function E(vel1,vel2,yy1,yy2) result(ene)
        implicit none
        real(8) :: vel1,vel2,yy1,yy2,ene
        ene = g*m1*yy1+g*m2*yy2+0.5d0*(m1*vel1**2.d0+m2*vel2**2.d0)
    end function E

end program double_pendulum




