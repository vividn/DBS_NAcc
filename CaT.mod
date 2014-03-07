TITLE calcium T channel for GPi neuron model

COMMENT

 Low threshold calcium channel (T-type), Wang et al. 1991 
 & Coulter et al 1989.  The original data was recorded at 22-24degC.
 Implementation from Gillies2006.

 Q10=1.52 --> rate_k=exp(log(Q10)*((1/296)-(1/309))/((1/292)-(1/302)))=1.69

ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaT
	USEION ca READ cai,cao,eca WRITE ica
	RANGE gcaT, iCaT
	GLOBAL rate_k,gmax_k
}

PARAMETER {
    v (mV)
	dt (ms)
	gcaT  = 0.001 (mho/cm2)
	iCaT  = 0.0 (mA/cm2)
	eca
	cai
	cao
	celsius
}

STATE {
    r s d
}

ASSIGNED { 
    ica (mA/cm2)
	ralpha (/ms)
	rbeta (/ms)
	salpha (/ms)
	sbeta (/ms)
	dalpha (/ms)
	dbeta (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica  = (gcaT*gmax_k)*r*r*r*s*ghkg(v,cai,cao,2)
	iCaT = ica
}

UNITSOFF

INITIAL {
    rate_k = 1.69
    gmax_k = 1.69
    settables(v)
    r = ralpha/(ralpha+rbeta)
    s = (salpha*(dbeta+dalpha) - (salpha*dbeta))/((salpha+sbeta)*(dalpha+dbeta) - (salpha*dbeta))
	d = (dbeta*(salpha+sbeta) - (salpha*dbeta))/((salpha+sbeta)*(dalpha+dbeta) - (salpha*dbeta))
}

DERIVATIVE states {  
	settables(v)      :Computes state variables at the current v and dt.
	r' = ((ralpha*(1-r)) - (rbeta*r))
	d' = ((dbeta*(1-s-d)) - (dalpha*d))
	s' = ((salpha*(1-s-d)) - (sbeta*s))
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
                          :Voltage shift (for temp effects) of -1.9278 added
    LOCAL   bd
    TABLE ralpha, rbeta, salpha, sbeta, dalpha, dbeta DEPEND celsius FROM -100 TO 100 WITH 400

	:"r" CaT activation system
	ralpha = rate_k * 1.0/(1.7+exp(-(v + 26.2722)/13.5))
	rbeta  = rate_k * exp(-(v + 61.0722)/7.8)/(exp(-(v + 26.8722)/13.1)+1.7)

    :"s" CaT fast inactivation system
    salpha = rate_k * exp(-(v + 158.3722)/17.8)
    sbeta  = rate_k * (sqrt(0.25+exp((v + 81.5722)/6.3))-0.5) * (exp(-(v + 158.3722)/17.8))

	:"d" CaT slow inactivation system
	bd     = sqrt(0.25+exp((v + 81.5722)/6.3))
    dalpha = rate_k * (1.0+exp((v + 35.4722)/30.0))/(240.0*(0.5+bd))
    dbeta  = rate_k * (bd-0.5)*dalpha
}

UNITSON

INCLUDE "ghk.inc"

