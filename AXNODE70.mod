TITLE Motor Axon Node channels
:
: I modified Cameron's AXNODE75.mod to have same parameters as axon_test_original.hoc (v_init=-70)
:
: Fast Na+, Persistant Na+, Slow K+, and Leakage currents 
: responsible for nodal action potential
: Iterative equations H-H notation rest = -75 mV
:


NEURON {
	SUFFIX axnode70
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT inap
	NONSPECIFIC_CURRENT ik
	NONSPECIFIC_CURRENT il
	RANGE gnapbar, gnabar, gkbar, gl, ena, ek, el
	RANGE mp_inf, m_inf, h_inf, s_inf
	RANGE tau_mp, tau_m, tau_h, tau_s
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	celsius		(degC)
	gnapbar = 0.05	(mho/cm2)
	gnabar	= 3.0	(mho/cm2)
	gkbar   = 0.07 	(mho/cm2)
	gl	= 0.005 (mho/cm2)
	ena     = 45.0  (mV)
	ek      = -95.0 (mV)
	el	= -60.0 (mV)
	vshift=10
	vtraub=-80
	ampA = 0.01
	ampB = 27
	ampC = 10.2
	bmpA = 0.00025
	bmpB = 34
	bmpC = 10
	amA = 1.86
	amB = 21.4
	amC = 10.3
	bmA = 0.086
	bmB = 25.7
	bmC = 9.16
	ahA = 0.062
	ahB = 114.0
	ahC = 11.0
	bhA = 2.3
	bhB = 31.8
	bhC = 13.4
	asA = 0.3
	asB = -27
	asC = -5
	bsA = 0.03
	bsB = 10
	bsC = -1
}

STATE {
	mp m h s
}

ASSIGNED {
	v       (mV)
	inap    (mA/cm2)
	ina	(mA/cm2)
	ik      (mA/cm2)
	il      (mA/cm2)
	mp_inf
	m_inf
	h_inf
	s_inf
	tau_mp
	tau_m
	tau_h
	tau_s
	q10_1
	q10_2
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	inap = gnapbar * mp*mp*mp * (v - ena)
	ina = gnabar * m*m*m*h * (v - ena)
	ik   = gkbar * s * (v - ek)
	il   = gl * (v - el)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	mp'= (mp_inf - mp) / tau_mp
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	s' = (s_inf - s) / tau_s
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:

	q10_1 = 2.2 ^ ((celsius-20)/ 10 )
	q10_2 = 2.9 ^ ((celsius-20)/ 10 )
	q10_3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	mp = mp_inf
	m = m_inf
	h = h_inf
	s = s_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vshift

	a = q10_1*vtrap1(v2)
	b = q10_1*vtrap2(v2)
	tau_mp = 1 / (a + b)
	mp_inf = a / (a + b)

	a = q10_1*vtrap6(v2)
	b = q10_1*vtrap7(v2)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	a = q10_2*vtrap8(v2)
	b = q10_2*bhA / (1 + Exp(-(v2+bhB)/bhC))
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)

	a = q10_3*asA / (Exp((v2-vtraub+asB)/asC) + 1) 
	b = q10_3*bsA / (Exp((v2-vtraub+bsB)/bsC) + 1)
	tau_s = 1 / (a + b)
	s_inf = a / (a + b)
}

FUNCTION vtrap(x) {
	if (x < -50) {
		vtrap = 0
	}else{
		vtrap = bsA / (Exp((x+bsB)/bsC) + 1)
	}
}

FUNCTION vtrap1(x) {
	if (fabs((x+ampB)/ampC) < 1e-6) {
		vtrap1 = ampA*ampC
	}else{
		vtrap1 = (ampA*(x+ampB)) / (1 - Exp(-(x+ampB)/ampC))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+bmpB)/bmpC) < 1e-6) {
		vtrap2 = -bmpA*bmpC
	}else{
		vtrap2 = (bmpA*(-(x+bmpB))) / (1 - Exp((x+bmpB)/bmpC))
	}
}

FUNCTION vtrap6(x) {
	if (fabs((x+amB)/amC) < 1e-6) {
		vtrap6 = amA*amC
	}else{
		vtrap6 = (amA*(x+amB)) / (1 - Exp(-(x+amB)/amC))
	}
}

FUNCTION vtrap7(x) {
	if (fabs((x+bmB)/bmC) < 1e-6) {
		vtrap7 = -bmA*bmC
	}else{
		vtrap7 = (bmA*(-(x+bmB))) / (1 - Exp((x+bmB)/bmC))
	}
}

FUNCTION vtrap8(x) {
	if (fabs((x+ahB)/ahC) < 1e-6) {
		vtrap8 = -ahA*ahC
	}else{
		vtrap8 = (ahA*(-(x+ahB))) / (1 - Exp((x+ahB)/ahC)) 
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON