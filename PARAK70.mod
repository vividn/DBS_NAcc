TITLE Paranode Axon channels
:
: I modified Cameron's PARAK75.mod to have same parameters as axon_test_original.hoc (v_init=-70)
:
: fast K+ currents 
: Iterative equations H-H notation rest = -75 mV
:


NEURON {
	SUFFIX parak70
	NONSPECIFIC_CURRENT ik
	RANGE gkbar, ek
	RANGE n_inf
	RANGE tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	celsius		(degC)
	gkbar   = 0.02	(mho/cm2)
	ek      = -95.0 (mV)
	vshift = 10	(mV)
	anA = 0.00798
	anB = 93.2
	anC = 1.1
	bnA = 0.0142
	bnB = 76
	bnC = 10.5
}

STATE {
	n
}

ASSIGNED {
	v               (mV)
	ik      (mA/cm2)
	n_inf
	tau_n
	q10
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik   = gkbar * n*n*n*n * (v - ek)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:

	q10 = 3.0 ^ ((celsius-20)/ 10 )

	evaluate_fct(v)
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v-vshift

	a = q10*vtrap1(v2)
	b = q10*vtrap2(v2)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)

}

FUNCTION vtrap1(x) {
	if (fabs((x+anB)/anC) < 1e-6) {
		vtrap1 = anA*anC
	}else{
		vtrap1 = (anA*(x+anB)) / (1 - Exp(-(x+anB)/anC))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+bnB)/bnC) < 1e-6) {
		vtrap2 = -bnA*bnC
	}else{
		vtrap2 = (bnA*(-(x+bnB))) / (1 - Exp((x+bnB)/bnC))
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