# Tables generated for v0.2.0

Unrounded source values are in the adjacent CSV files.

## Modal verification

| model | mode | natural frequency hz | damped frequency hz | damping ratio | reference | equation error |
|:--|:--|:--|:--|:--|:--|:--|
| Planar torsional pendulum | 1 | 0.7750535255835627 | 0.7735814377863248 | 0.061604110363369734 | closed form | 2.303098026953073e-16 |
| Spatial torsional pendulum | 1 | 0.7796968012336759 | 0.7775887701761225 | 0.07348469228349536 | closed form | 1.0925618424287063e-16 |
| Spatial three-link pendulum | 1 | 0.37018333469035547 | 0.3701819300577104 | 0.0027547838801786156 | independent state matrix | 1.23404758491105e-15 |
| Spatial three-link pendulum | 2 | 1.125806559191502 | 1.125344627562513 | 0.028643587243083975 | independent state matrix | 8.938355511080724e-16 |
| Spatial three-link pendulum | 3 | 2.4213613788934576 | 2.412938174825621 | 0.08333852870440284 | independent state matrix | 9.583726546659211e-16 |

## Open pendulum chains

| size | variables | states | jacobian nonzeros | load seconds | run seconds | allocated mib | accepted steps | rejected steps | corrector failures | residual evaluations | jacobian evaluations | numerical factorizations | symbolic factorizations | newton iterations | maximum order | maximum joint error |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 10 | 110 | 10 | 498 | 0.001491458 | 0.010892 | 49.409027099609375 | 95 | 2 | 0 | 327 | 20 | 20 | 1 | 229 | 5 | 3.7317369404714316e-11 |
| 25 | 275 | 25 | 1278 | 0.002555833 | 0.031159792 | 148.9141387939453 | 123 | 2 | 0 | 392 | 25 | 25 | 1 | 266 | 5 | 8.351741858013829e-12 |
| 50 | 550 | 50 | 2578 | 0.006186417 | 0.102637958 | 349.728271484375 | 147 | 2 | 0 | 454 | 30 | 30 | 1 | 304 | 5 | 1.7682107813646217e-12 |

## State choices

| state choice | run seconds | allocated mib | variables | selected states | accepted steps | rejected steps | residual evaluations | numerical factorizations | symbolic factorizations | newton iterations | corrector failures |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| Automatic QR | 0.922758 | 3304.1009063720703 | 120 | 10 | 6171 | 75 | 22395 | 1312 | 68 | 16148 | 0 |
| Preferred body angular velocities | 0.943460209 | 3395.6956329345703 | 120 | 10 | 6254 | 135 | 22960 | 1353 | 75 | 16570 | 0 |
| Preferred relative joint angular velocities | 1.126452375 | 3961.7739868164062 | 150 | 10 | 7130 | 99 | 24139 | 1516 | 76 | 16909 | 0 |

## Closed parallelogram chains

| size | variables | states | jacobian nonzeros | load seconds | run seconds | allocated mib | accepted steps | rejected steps | corrector failures | residual evaluations | jacobian evaluations | numerical factorizations | symbolic factorizations | newton iterations | maximum order | maximum joint error |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 10 | 262 | 1 | 1254 | 0.003001833 | 0.03227225 | 132.2600555419922 | 69 | 1 | 0 | 312 | 23 | 23 | 3 | 241 | 5 | 3.3671063623353387e-10 |
| 25 | 637 | 1 | 3084 | 0.009355209 | 0.120409 | 347.8019561767578 | 69 | 1 | 0 | 328 | 23 | 23 | 3 | 257 | 5 | 3.450734786612142e-10 |
| 50 | 1262 | 1 | 6134 | 0.037339208 | 0.286355125 | 913.0090942382812 | 69 | 1 | 0 | 332 | 23 | 23 | 2 | 261 | 5 | 3.2439614617969353e-10 |

## Rotor trains

| rotors | variables | states | run seconds | allocated mib | accepted steps | rejected steps | corrector failures | numerical factorizations | symbolic factorizations | minimum frequency rad s | maximum frequency rad s | maximum angle error rad | maximum angular velocity error rad s | final energy ratio |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 10 | 120 | 10 | 0.066172083 | 294.36912536621094 | 908 | 11 | 0 | 185 | 1 | 9.45269221962589 | 125.07830525830259 | 5.960084810303912e-7 | 6.126580618870392e-5 | 0.03331871073611405 |
| 25 | 300 | 25 | 0.159465917 | 692.2382202148438 | 867 | 1 | 0 | 175 | 1 | 3.8953010286295986 | 126.25119436019402 | 5.93470456656392e-7 | 6.697208880357408e-5 | 0.03355567287771608 |
| 50 | 600 | 50 | 0.315491667 | 1286.6349029541016 | 827 | 0 | 0 | 166 | 1 | 1.9671658964347356 | 126.42992041866022 | 7.76525372857273e-7 | 8.753199085376195e-5 | 0.03351595300267536 |

## Rotating flexible blades

| segments | target speed rad s | actual speed rad s | static tip height m | final tip height m | first frequency hz | maximum modal equation error | accepted steps | rejected steps | corrector failures | state reselections |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 4 | 0.0 | 0.0 | -0.46555182337184176 | -0.46555182337183965 | 0.9108915437801168 | 8.394097761562787e-10 | 616 | 0 | 0 | 0 |
| 4 | 7.5 | 7.5 | -0.46555182337183965 | -0.1529916699413396 | 1.5639807403720989 | 6.718931308498814e-12 | 699 | 4 | 0 | 0 |
| 8 | 7.5 | 7.5 | -0.4655967016235369 | -0.15108658822368312 | 1.5746731911605498 | 1.2607173458315183e-8 | 710 | 24 | 0 | 19 |

## Bouncing-ball soft restarts

| balls | policy | variables | run seconds | allocated mib | accepted steps | rejected steps | corrector failures | root evaluations | events | history restarts | maximum sampled penetration m |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 10 | soft | 140 | 0.263747959 | 923.6025238037109 | 1578 | 220 | 0 | 3027 | 45 | 45 | 0.05642228016700793 |
| 25 | soft | 350 | 1.100556084 | 4548.424880981445 | 2442 | 493 | 0 | 6141 | 112 | 112 | 0.05642163244266403 |
| 50 | soft | 700 | 3.4798255 | 16235.706848144531 | 3471 | 874 | 0 | 10841 | 224 | 224 | 0.05642188818590768 |

## Large Van implementation microbenchmark

| implementation | active variables | run seconds | allocated mib | allocations | gc seconds | static iterations | accepted steps | rejected steps | residual evaluations | jacobian evaluations | numerical factorizations | symbolic factorizations | newton iterations | final maximum |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 0.2.0 current implementation | 1319 | 0.59309325 | 902.0372009277344 | 10191316 | 0.029083041 | 28 | 57 | 0 | 524 | 50 | 50 | 1 | 466 | 11961.755354970339 |

## Large Van 10-second maneuver

| maneuver | simulated seconds | output samples | run seconds | simulated seconds per run second | allocated mib | gc seconds | static iterations | initial ground speed m s | minimum ground speed m s | minimum speed time s | final ground speed m s | initial body forward speed m s | final body forward speed m s | accepted steps | rejected steps | residual evaluations | jacobian evaluations | numerical factorizations | symbolic factorizations | newton iterations | corrector failures | state reselections |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 180 degree steering-wheel sweep | 10.0 | 601 | 11.099682625 | 0.9009266605044034 | 18740.827407836914 | 1.379377031 | 28 | 30.0 | 1.276030550103443 | 5.516666666666667 | 1.601331675546403 | 29.99672582444076 | 1.5943917975661477 | 1292 | 9 | 14860 | 1248 | 1248 | 6 | 13565 | 7 | 2 |
