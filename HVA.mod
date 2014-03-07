TITLE calcium HVA channels for GPi neuron model

COMMENT

 High threshold calcium channel (N/L-type), Brown et al. 1993. & Fox
 et al. (1989).  Both done at temperature 22degC.  Implemented from
 Gillies2006.

 Adding CaL [Ca]i dependent inactivation.  This is only for the L-type
 component, and is called inactivation variable 'h'.  

 Q10=1.95 --> rate_k=exp(log(Q10)*((1/295)-(1/309))/((1/293)-(1/303)))=2.49

ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mV) = (millivolt)
	(mA) = (milliamp)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX HVA
	USEION ca READ cai,cao,eca WRITE ica
	RANGE gcaN, gcaL, iNCa, iLCa
	GLOBAL inactLtau,inactLmax,rate_k,gmax_k
}

PARAMETER {
    v (mV)
	dt (ms)
    gcaL  = 0.002 (mho/cm2)
	gcaN  = 0.012 (mho/cm2)
	iNCa  = 0.0 (mA/cm2)
	iLCa  = 0.0 (mA/cm2)
	inactLtau = 1220.0 (ms)
	inactLmax = 0.529
	eca
	cai
	cao
	celsius
}

STATE {
    q u h
}

ASSIGNED { 
    ica (mA/cm2)
	qinf
	uinf
	hinf
	qtau (ms)
	utau (ms)
	htau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	LOCAL vghk
	SOLVE states METHOD cnexp
	vghk = ghkg(v,cai,cao,2)
	iNCa = gmax_k*(gcaN * u)*q*q*vghk
	iLCa = gmax_k*(gcaL)*q*q*h*vghk
	ica  = iNCa + iLCa
}

INITIAL {
    rate_k = 2.49
	gmax_k = 2.49
    settables(v)
	q = qinf
	u = uinf
	setCadepLinact(cai)
	h = hinf
}

DERIVATIVE states {  
	settables(v)  
	q' = (qinf-q)/qtau
	u' = (uinf-u)/utau
	setCadepLinact(cai)
	h' = (hinf-h)/htau
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
                          :Voltage shifts (for temp effects) of -8.25 and -14.67 added respt.
    TABLE qinf, qtau, uinf, utau DEPEND celsius FROM -100 TO 100 WITH 400

    :"q" N/L Ca activation system
    qinf   = 1.0/(1.0 + exp((-16.3547869 - v)/11.3))
    qtau   = (1.25/(cosh(-0.031 * (v + 28.8547869)))) /rate_k

    :"u" N inactivation system - voltage dependent.
    uinf   = 1.0/(1.0 + exp((v + 45.3326653)/12.5))
    utau   = (98.0 + cosh(0.021*(24.7673347-v))) /rate_k
}

PROCEDURE setCadepLinact(cai) { : set Ca dependent L-type calcium channel inactivation
    :"h" L inactivation system - [Ca]i dependent.
	hinf   = inactLmax+((1.0-inactLmax)/(1.0 + exp((cai-0.7)/0.15)))
    htau   = inactLtau /rate_k
}

INCLUDE "ghk.inc"


