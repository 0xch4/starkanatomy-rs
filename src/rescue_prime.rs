//! Rescue Prime hashing.

use crate::{algebra::FieldElement, multivariate::MultiPolynomial, univariate::Polynomial};

pub const M: usize = 2;
pub const N: usize = 27;
const ALPHA: u128 = 3;
const ALPHAINV: u128 = 180331931428153586757283157844700080811;
const MDS: [[FieldElement; 2]; 2] = [
    [
        FieldElement::new(270497897142230380135924736767050121214),
        FieldElement::new(4),
    ],
    [
        FieldElement::new(270497897142230380135924736767050121205),
        FieldElement::new(13),
    ],
];
const MDSINV: [[FieldElement; 2]; 2] = [
    [
        FieldElement::new(210387253332845851216830350818816760948),
        FieldElement::new(60110643809384528919094385948233360270),
    ],
    [
        FieldElement::new(90165965714076793378641578922350040407),
        FieldElement::new(180331931428153586757283157844700080811),
    ],
];
const ROUND_CONSTANTS: [FieldElement; 108] = [
    FieldElement::new(174420698556543096520990950387834928928),
    FieldElement::new(109797589356993153279775383318666383471),
    FieldElement::new(228209559001143551442223248324541026000),
    FieldElement::new(268065703411175077628483247596226793933),
    FieldElement::new(250145786294793103303712876509736552288),
    FieldElement::new(154077925986488943960463842753819802236),
    FieldElement::new(204351119916823989032262966063401835731),
    FieldElement::new(57645879694647124999765652767459586992),
    FieldElement::new(102595110702094480597072290517349480965),
    FieldElement::new(8547439040206095323896524760274454544),
    FieldElement::new(50572190394727023982626065566525285390),
    FieldElement::new(87212354645973284136664042673979287772),
    FieldElement::new(64194686442324278631544434661927384193),
    FieldElement::new(23568247650578792137833165499572533289),
    FieldElement::new(264007385962234849237916966106429729444),
    FieldElement::new(227358300354534643391164539784212796168),
    FieldElement::new(179708233992972292788270914486717436725),
    FieldElement::new(102544935062767739638603684272741145148),
    FieldElement::new(65916940568893052493361867756647855734),
    FieldElement::new(144640159807528060664543800548526463356),
    FieldElement::new(58854991566939066418297427463486407598),
    FieldElement::new(144030533171309201969715569323510469388),
    FieldElement::new(264508722432906572066373216583268225708),
    FieldElement::new(22822825100935314666408731317941213728),
    FieldElement::new(33847779135505989201180138242500409760),
    FieldElement::new(146019284593100673590036640208621384175),
    FieldElement::new(51518045467620803302456472369449375741),
    FieldElement::new(73980612169525564135758195254813968438),
    FieldElement::new(31385101081646507577789564023348734881),
    FieldElement::new(270440021758749482599657914695597186347),
    FieldElement::new(185230877992845332344172234234093900282),
    FieldElement::new(210581925261995303483700331833844461519),
    FieldElement::new(233206235520000865382510460029939548462),
    FieldElement::new(178264060478215643105832556466392228683),
    FieldElement::new(69838834175855952450551936238929375468),
    FieldElement::new(75130152423898813192534713014890860884),
    FieldElement::new(59548275327570508231574439445023390415),
    FieldElement::new(43940979610564284967906719248029560342),
    FieldElement::new(95698099945510403318638730212513975543),
    FieldElement::new(77477281413246683919638580088082585351),
    FieldElement::new(206782304337497407273753387483545866988),
    FieldElement::new(141354674678885463410629926929791411677),
    FieldElement::new(19199940390616847185791261689448703536),
    FieldElement::new(177613618019817222931832611307175416361),
    FieldElement::new(267907751104005095811361156810067173120),
    FieldElement::new(33296937002574626161968730356414562829),
    FieldElement::new(63869971087730263431297345514089710163),
    FieldElement::new(200481282361858638356211874793723910968),
    FieldElement::new(69328322389827264175963301685224506573),
    FieldElement::new(239701591437699235962505536113880102063),
    FieldElement::new(17960711445525398132996203513667829940),
    FieldElement::new(219475635972825920849300179026969104558),
    FieldElement::new(230038611061931950901316413728344422823),
    FieldElement::new(149446814906994196814403811767389273580),
    FieldElement::new(25535582028106779796087284957910475912),
    FieldElement::new(93289417880348777872263904150910422367),
    FieldElement::new(4779480286211196984451238384230810357),
    FieldElement::new(208762241641328369347598009494500117007),
    FieldElement::new(34228805619823025763071411313049761059),
    FieldElement::new(158261639460060679368122984607245246072),
    FieldElement::new(65048656051037025727800046057154042857),
    FieldElement::new(134082885477766198947293095565706395050),
    FieldElement::new(23967684755547703714152865513907888630),
    FieldElement::new(8509910504689758897218307536423349149),
    FieldElement::new(232305018091414643115319608123377855094),
    FieldElement::new(170072389454430682177687789261779760420),
    FieldElement::new(62135161769871915508973643543011377095),
    FieldElement::new(15206455074148527786017895403501783555),
    FieldElement::new(201789266626211748844060539344508876901),
    FieldElement::new(179184798347291033565902633932801007181),
    FieldElement::new(9615415305648972863990712807943643216),
    FieldElement::new(95833504353120759807903032286346974132),
    FieldElement::new(181975981662825791627439958531194157276),
    FieldElement::new(267590267548392311337348990085222348350),
    FieldElement::new(49899900194200760923895805362651210299),
    FieldElement::new(89154519171560176870922732825690870368),
    FieldElement::new(265649728290587561988835145059696796797),
    FieldElement::new(140583850659111280842212115981043548773),
    FieldElement::new(266613908274746297875734026718148328473),
    FieldElement::new(236645120614796645424209995934912005038),
    FieldElement::new(265994065390091692951198742962775551587),
    FieldElement::new(59082836245981276360468435361137847418),
    FieldElement::new(26520064393601763202002257967586372271),
    FieldElement::new(108781692876845940775123575518154991932),
    FieldElement::new(138658034947980464912436420092172339656),
    FieldElement::new(45127926643030464660360100330441456786),
    FieldElement::new(210648707238405606524318597107528368459),
    FieldElement::new(42375307814689058540930810881506327698),
    FieldElement::new(237653383836912953043082350232373669114),
    FieldElement::new(236638771475482562810484106048928039069),
    FieldElement::new(168366677297979943348866069441526047857),
    FieldElement::new(195301262267610361172900534545341678525),
    FieldElement::new(2123819604855435621395010720102555908),
    FieldElement::new(96986567016099155020743003059932893278),
    FieldElement::new(248057324456138589201107100302767574618),
    FieldElement::new(198550227406618432920989444844179399959),
    FieldElement::new(177812676254201468976352471992022853250),
    FieldElement::new(211374136170376198628213577084029234846),
    FieldElement::new(105785712445518775732830634260671010540),
    FieldElement::new(122179368175793934687780753063673096166),
    FieldElement::new(126848216361173160497844444214866193172),
    FieldElement::new(22264167580742653700039698161547403113),
    FieldElement::new(234275908658634858929918842923795514466),
    FieldElement::new(189409811294589697028796856023159619258),
    FieldElement::new(75017033107075630953974011872571911999),
    FieldElement::new(144945344860351075586575129489570116296),
    FieldElement::new(261991152616933455169437121254310265934),
    FieldElement::new(18450316039330448878816627264054416127),
];

/// Rescue Prime hasher.
pub struct RescuePrime {}

impl RescuePrime {
    pub fn new() -> Self {
        Self {}
    }

    pub fn hash(&self, input: FieldElement) -> FieldElement {
        let mut state = [[input], [FieldElement::ZERO; M - 1]].concat();

        for r in 0..N {
            state = state.iter().map(|fe| fe.pow(ALPHA)).collect();

            let tmp: Vec<FieldElement> = MDS
                .iter()
                .map(|mds_row| {
                    mds_row
                        .iter()
                        .zip(state.iter())
                        .map(|(&mds, &fe)| mds * fe)
                        .sum()
                })
                .collect();
            state = tmp
                .iter()
                .enumerate()
                .map(|(i, &tmp_fe)| tmp_fe + ROUND_CONSTANTS[2 * r * M + i])
                .collect();

            state = state.iter().map(|fe| fe.pow(ALPHAINV)).collect();

            let tmp: Vec<FieldElement> = MDS
                .iter()
                .map(|mds_row| {
                    mds_row
                        .iter()
                        .zip(state.iter())
                        .map(|(&mds, &fe)| mds * fe)
                        .sum()
                })
                .collect();
            state = tmp
                .iter()
                .enumerate()
                .map(|(i, &tmp_fe)| tmp_fe + ROUND_CONSTANTS[2 * r * M + M + i])
                .collect();
        }

        state[0]
    }

    pub fn trace(&self, input: FieldElement) -> Vec<Vec<FieldElement>> {
        let mut trace = Vec::new();

        let mut state = [[input], [FieldElement::ZERO; M - 1]].concat();
        trace.push(state.clone());

        for r in 0..N {
            state = state.iter().map(|fe| fe.pow(ALPHA)).collect();

            let tmp: Vec<FieldElement> = MDS
                .iter()
                .map(|mds_row| {
                    mds_row
                        .iter()
                        .zip(state.iter())
                        .map(|(&mds, &fe)| mds * fe)
                        .sum()
                })
                .collect();
            state = tmp
                .iter()
                .enumerate()
                .map(|(i, &tmp_fe)| tmp_fe + ROUND_CONSTANTS[2 * r * M + i])
                .collect();

            state = state.iter().map(|fe| fe.pow(ALPHAINV)).collect();

            let tmp: Vec<FieldElement> = MDS
                .iter()
                .map(|mds_row| {
                    mds_row
                        .iter()
                        .zip(state.iter())
                        .map(|(&mds, &fe)| mds * fe)
                        .sum()
                })
                .collect();
            state = tmp
                .iter()
                .enumerate()
                .map(|(i, &tmp_fe)| tmp_fe + ROUND_CONSTANTS[2 * r * M + M + i])
                .collect();

            trace.push(state.clone());
        }

        trace
    }

    fn round_constants_polynomials(
        &self,
        omicron: FieldElement,
    ) -> (Vec<MultiPolynomial>, Vec<MultiPolynomial>) {
        let mut first_step_constants = Vec::new();
        for i in 0..M {
            let domain = (0..N)
                .map(|r| (omicron.pow(r as _), ROUND_CONSTANTS[2 * r * M + i]))
                .collect();
            first_step_constants.push(MultiPolynomial::lift_univariate(
                &Polynomial::interpolate(&domain),
                1,
                0,
            ));
        }

        let mut second_step_constants = Vec::new();
        for i in 0..M {
            let domain = (0..N)
                .map(|r| (omicron.pow(r as _), ROUND_CONSTANTS[2 * r * M + M + i]))
                .collect();
            second_step_constants.push(MultiPolynomial::lift_univariate(
                &Polynomial::interpolate(&domain),
                1,
                0,
            ));
        }

        (first_step_constants, second_step_constants)
    }

    pub fn boundary_constraints(&self, output: FieldElement) -> Vec<(usize, usize, FieldElement)> {
        vec![(0, 1, FieldElement::ZERO), (N, 0, output)]
    }

    pub fn transition_constraints(&self, omicron: FieldElement) -> Vec<MultiPolynomial> {
        let (first_step_constants, second_step_constants) =
            self.round_constants_polynomials(omicron);

        let variables = MultiPolynomial::univariate_variables(1 + 2 * M);
        let previous_state = variables[1..(1 + M)].to_vec();
        let next_state = variables[(1 + M)..(1 + 2 * M)].to_vec();

        let mut air = Vec::new();
        for i in 0..M {
            let mut lhs = MultiPolynomial::ZERO;
            for k in 0..M {
                lhs = lhs.add(
                    &MultiPolynomial::constant(MDS[i][k], 1).mul(&previous_state[k].pow(ALPHA)),
                );
            }
            lhs = lhs.add(&first_step_constants[i]);

            let mut rhs = MultiPolynomial::ZERO;
            for k in 0..M {
                rhs = rhs.add(
                    &MultiPolynomial::constant(MDSINV[i][k], 1)
                        .mul(&next_state[k].sub(&second_step_constants[k])),
                );
            }
            rhs = rhs.pow(ALPHA);

            air.push(lhs.sub(&rhs));
        }

        air
    }
}

#[cfg(test)]
mod tests {
    use rand::{rngs::OsRng, Rng};

    use super::*;

    #[test]
    fn hash() {
        let rp = RescuePrime::new();
        assert_eq!(
            rp.hash(FieldElement::ONE),
            FieldElement::new(244180265933090377212304188905974087294)
        );
        assert_eq!(
            rp.hash(FieldElement::new(57322816861100832358702415967512842988)),
            FieldElement::new(89633745865384635541695204788332415101)
        );
    }

    #[test]
    fn trace() {
        let rp = RescuePrime::new();
        let a = FieldElement::new(57322816861100832358702415967512842988);
        let b = FieldElement::new(89633745865384635541695204788332415101);
        let trace = rp.trace(a);
        assert_eq!(trace[0][0], a);
        assert_eq!(trace[trace.len() - 1][0], b);
    }

    #[test]
    fn constraints() {
        let rp = RescuePrime::new();
        let input = FieldElement::new(57322816861100832358702415967512842988);
        let output = rp.hash(input);
        let trace = rp.trace(input);

        for &condition in rp.boundary_constraints(output).iter() {
            let (cycle, element, value) = condition;
            assert_eq!(trace[cycle][element], value);
        }

        let omicron = FieldElement::primitive_nth_root(1 << 119);
        let transition_constraints = rp.transition_constraints(omicron);
        for o in 0..(trace.len() - 1) {
            for air_poly in transition_constraints.iter() {
                let previous_state = vec![trace[o][0], trace[o][1]];
                let next_state = vec![trace[o + 1][0], trace[o + 1][1]];
                let point = [vec![omicron.pow(o as _)], previous_state, next_state].concat();
                assert_eq!(air_poly.evaluate(&point), FieldElement::ZERO);
            }
        }
    }

    #[test]
    fn invalid() {
        let rp = RescuePrime::new();
        let input = FieldElement::new(57322816861100832358702415967512842988);
        let output = rp.hash(input);
        let mut trace = rp.trace(input);
        let omicron = FieldElement::primitive_nth_root(1 << 119);
        let transition_constraints = rp.transition_constraints(omicron);
        let boundary_constraints = rp.boundary_constraints(output);

        for trial in 0..10 {
            let register_index = if trial > 0 { OsRng.gen_range(0..M) } else { 1 };
            let cycle_index = if trial > 0 {
                OsRng.gen_range(0..(N + 1))
            } else {
                22
            };
            let mut arr = [0u8; 16];
            OsRng.fill(&mut arr);
            let value = if trial > 0 {
                FieldElement::sample(arr)
            } else {
                FieldElement::new(17274817952119230544216945715808633996)
            };
            if value == FieldElement::ZERO {
                continue;
            }

            trace[cycle_index][register_index] = trace[cycle_index][register_index] + value;
            let mut error_got_noticed = false;

            for &condition in boundary_constraints.iter() {
                if error_got_noticed {
                    break;
                }
                let (cycle, element, value) = condition;
                if trace[cycle][element] != value {
                    error_got_noticed = true;
                }
            }

            for o in 0..trace.len() - 1 {
                if error_got_noticed {
                    break;
                }
                for air_poly in transition_constraints.iter() {
                    let previous_state = vec![trace[o][0], trace[o][1]];
                    let next_state = vec![trace[o + 1][0], trace[o + 1][1]];
                    let point = [vec![omicron.pow(o as _)], previous_state, next_state].concat();
                    if air_poly.evaluate(&point) != FieldElement::ZERO {
                        error_got_noticed = true;
                    }
                }
            }

            assert!(error_got_noticed);
            trace[cycle_index][register_index] = trace[cycle_index][register_index] - value;
        }
    }
}
