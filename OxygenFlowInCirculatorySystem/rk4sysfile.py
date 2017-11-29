from numpy import diff,any,arange,append
from sys import exit


def rk4sys(dydt,tspan,y0,h,*args):
    '''
    rk4sys: 4th-order Runge-Kutta for a system of ODEs
    
    [t,y] = rk4sys(dydt,tspan,y0,h,p1,p2,...): integrates
    a system of ODEs with 4th-order RK method
    
    input:
    dydt = name of the Python file that evaluates the ODEs
    tspan = [ti,tf]; initial and final times with output
    generated at interval of h, or
    = [t0 t1 ... tf]; specific times where solution output
    y0 = initial values of dependent variables
    h = step size
    p1,p2,... = additional parameters used by dydt
    
    output:
    tp = vector of independent variable
    yp = vector of solution for dependent variables
    '''
    if any(diff(tspan)) <= 0.:
        print '\ntspan not ascending order.\n'
        exit(0)
    n = len(tspan)
    ti = tspan[0]
    tf = tspan[n-1]
    if n == 2:
        t = arange(ti,tf+h,h)
        n = len(t)
        if t[n-1] < tf:
            t = append(t,tf)
            n = n+1
    else:
        t = tspan
    tt = ti
    y = y0
    np = 0
    tp = tt
    yp = y
    i = 0
    timekeep=0
    while(1):
        tend = t[np+1]
        hh = t[np+1]-t[np]
        if hh > h:
            hh = h
        while(1):
            if tt+hh > tend:
                hh = tend-tt
            k1 = dydt(tt,y[i],args)
            ymid = y[i]+k1*hh/2.
            k2 = dydt(tt+hh/2.,ymid[0],args)
            ymid = y[i]+k2*hh/2
            k3 = dydt(tt+hh/2.,ymid[0],args)
            yend = y[i]+k3*hh
            k4 = dydt(tt+hh,yend[0],args)
            phi = (k1+2.*(k2+k3)+k4)/6.
            y = append(y,y[i]+phi*hh,axis=0)
            tt = tt+hh
            i = i+1
            timekeep=timekeep+1
            if timekeep%10000 == 0:
                print 'Percent Completed:',tt/tf*100
            if tt >= tend:
                break
        np = np+1
        tp = append(tp,tt)
        yp = append(yp,[y[i]],axis=0)
        if tt >= tf:
            break
    return [tp,yp]
