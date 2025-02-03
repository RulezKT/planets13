package planets13

import (
	"fmt"
	"math"
	"path/filepath"

	"github.com/RulezKT/floatsfile"
	"github.com/RulezKT/types"
)

const (
	MOON_FILE       = "moon.bin"
	SUN_FILE        = "sun.bin"
	MER_FILE        = "mercury.bin"
	VEN_FILE        = "venus.bin"
	MAR_FILE        = "mars.bin"
	JUP_FILE        = "jupyter.bin"
	SAT_FILE        = "saturn.bin"
	URA_FILE        = "uranus.bin"
	NEP_FILE        = "neptune.bin"
	PLU_FILE        = "pluto.bin"
	EARTH_FILE      = "earth.bin"
	EARTH_BARY_FILE = "earth_bary.bin"
	CERES_FILE      = "ceres.bin"
	CHIRON_FILE     = "chiron.bin"

	EARTH_LENGTH      = 752801
	MOON_LENGTH       = 752801
	SUN_LENGTH        = 160685
	MER_LENGTH        = 403788
	VEN_LENGTH        = 146848
	MAR_LENGTH        = 80325
	JUP_LENGTH        = 59670
	SAT_LENGTH        = 52785
	URA_LENGTH        = 45900
	NEP_LENGTH        = 45900
	PLU_LENGTH        = 45900
	EARTH_BARY_LENGTH = 188149
	CERES_LENGTH      = 457984
	CHIRON_LENGTH     = 458880
)

type Pl13 struct {
	moon   []float64
	sun    []float64
	mer    []float64
	ven    []float64
	mar    []float64
	jup    []float64
	sat    []float64
	ura    []float64
	nep    []float64
	plu    []float64
	earth  []float64
	earthb []float64
	chiron []float64
	ceres  []float64
}

func (pl *Pl13) Load(dir string) {

	pl.earth = floatsfile.LoadBinary(filepath.Join(dir, EARTH_FILE), EARTH_LENGTH)
	pl.moon = floatsfile.LoadBinary(filepath.Join(dir, MOON_FILE), MOON_LENGTH)
	pl.sun = floatsfile.LoadBinary(filepath.Join(dir, SUN_FILE), SUN_LENGTH)
	pl.mer = floatsfile.LoadBinary(filepath.Join(dir, MER_FILE), MER_LENGTH)
	pl.ven = floatsfile.LoadBinary(filepath.Join(dir, VEN_FILE), VEN_LENGTH)
	pl.mar = floatsfile.LoadBinary(filepath.Join(dir, MAR_FILE), MAR_LENGTH)
	pl.jup = floatsfile.LoadBinary(filepath.Join(dir, JUP_FILE), JUP_LENGTH)
	pl.sat = floatsfile.LoadBinary(filepath.Join(dir, SAT_FILE), SAT_LENGTH)
	pl.ura = floatsfile.LoadBinary(filepath.Join(dir, URA_FILE), URA_LENGTH)
	pl.nep = floatsfile.LoadBinary(filepath.Join(dir, NEP_FILE), NEP_LENGTH)
	pl.plu = floatsfile.LoadBinary(filepath.Join(dir, PLU_FILE), PLU_LENGTH)
	pl.earthb = floatsfile.LoadBinary(filepath.Join(dir, EARTH_BARY_FILE), EARTH_BARY_LENGTH)
	pl.ceres = floatsfile.LoadBinary(filepath.Join(dir, CERES_FILE), CERES_LENGTH)
	pl.chiron = floatsfile.LoadBinary(filepath.Join(dir, CHIRON_FILE), CHIRON_LENGTH)
}

func (pl *Pl13) Calc(seconds float64) []types.Position {

	// earth barycenter
	earthBary := calcPos11(seconds, -3_156_062_400, 41, 1_382_400, 12, pl.earthb)

	// earth
	earthMoon := calcPos11(seconds, -3_157_963_200, 41, 345_600, 12, pl.earth)

	// sun
	sun := calcPos11(seconds, -3_157_444_800, 35, 1_382_400, 10, pl.sun)

	// chiron
	chiron := calcChiron(seconds, pl.chiron)
	chiron = chiron.Minus(&sun)
	chiron = chiron.Minus(&earthMoon)
	chiron = chiron.Minus(&earthBary)

	// ceres
	ceres := calcCeres(seconds, pl.ceres)
	ceres = ceres.Minus(&sun)
	ceres = ceres.Minus(&earthMoon)
	ceres = ceres.Minus(&earthBary)

	// sun
	sun = sun.Minus(&earthMoon)
	sun = sun.Minus(&earthBary)

	// moon
	moon := calcPos11(seconds, -3_157_963_200, 41, 345_600, 12, pl.moon)
	moon = moon.Minus(&earthBary)

	// mercury
	merc := calcPos11(seconds, -3_155_716_800, 44, 691_200, 13, pl.mer)
	merc = merc.Minus(&earthMoon)
	merc = merc.Minus(&earthBary)

	// venus
	ven := calcPos11(seconds, -3_156_062_400, 32, 1_382_400, 9, pl.ven)
	ven = ven.Minus(&earthMoon)
	ven = ven.Minus(&earthBary)

	// mars
	mars := calcPos11(seconds, -3_156_753_600, 35, 2_764_800, 10, pl.mar)
	mars = mars.Minus(&earthMoon)
	mars = mars.Minus(&earthBary)

	// jupiter
	jup := calcPos11(seconds, -3_156_753_600, 26, 2_764_800, 7, pl.jup)
	jup = jup.Minus(&earthMoon)
	jup = jup.Minus(&earthBary)

	// saturn
	sat := calcPos11(seconds, -3_156_753_600, 23, 2_764_800, 6, pl.sat)
	sat = sat.Minus(&earthMoon)
	sat = sat.Minus(&earthBary)

	// uranus
	uran := calcPos11(seconds, -3_156_753_600, 20, 2_764_800, 5, pl.ura)
	uran = uran.Minus(&earthMoon)
	uran = uran.Minus(&earthBary)

	// neptune
	nep := calcPos11(seconds, -3_156_753_600, 20, 2_764_800, 5, pl.nep)
	nep = nep.Minus(&earthMoon)
	nep = nep.Minus(&earthBary)

	// pluto
	plu := calcPos11(seconds, -3_156_753_600, 20, 2_764_800, 5, pl.plu)
	plu = plu.Minus(&earthMoon)
	plu = plu.Minus(&earthBary)

	return []types.Position{
		// earth
		{X: 0, Y: 0, Z: 0},

		// moon
		moon,

		// sun
		sun,

		// mercury
		merc,

		// venus
		ven,

		// mars
		mars,

		// ceres
		ceres,

		// jupiter
		jup,

		// saturn
		sat,

		// chiron
		chiron,

		// uranus
		uran,

		// neptune
		nep,

		// pluto
		plu,
	}

}

func calcPos11(seconds float64, startTime int64, rsize int, intlen int, order int, arrPtr []float64) types.Position {

	deg := order + 1

	offset := math.Floor((seconds-float64(startTime))/float64(intlen)) * float64(rsize)
	data := arrPtr[int64(offset) : int64(offset)+int64(rsize)]
	tau := (seconds - data[0]) / data[1]

	return types.Position{
		X: chebyshev(order, tau, data[2:2+deg]),
		Y: chebyshev(order, tau, data[2+deg:2+2*deg]),
		Z: chebyshev(order, tau, data[2+2*deg:2+3*deg]),
	}
}

func (pl *Pl13) CalcSun(seconds float64) types.Position {
	earthBary := calcPos11(seconds, -3_156_062_400, 41, 1_382_400, 12, pl.earthb)
	earthMoon := calcPos11(seconds, -3_157_963_200, 41, 345_600, 12, pl.earth)

	sun := calcPos11(seconds, -3_157_444_800, 35, 1_382_400, 10, pl.sun)
	sun = sun.Minus(&earthMoon)
	sun = sun.Minus(&earthBary)

	return sun

}

func chebyshev(order int, x float64, data []float64) float64 {

	// Evaluate a Chebyshev polynomial
	var bk float64
	two_x := 2 * x
	bkp2 := data[order]
	bkp1 := two_x*bkp2 + data[order-1]

	for n := order - 2; n > 0; n-- {
		bk = data[n] + two_x*bkp1 - bkp2
		bkp2 = bkp1
		bkp1 = bk
	}
	return data[0] + x*bkp1 - bkp2
}

type FileRecords struct {
	rec_start_addr int
	seg_start_time float64
	seg_last_time  float64
	int_len        float64
	rec_last_addr  int
}

func calcChiron(seconds float64, chiron []float64) types.Position {

	var CHIRON_FILE_RECORDS = []FileRecords{
		{
			// SUN
			rec_start_addr: 0,
			seg_start_time: 0,
			seg_last_time:  0,
			int_len:        0,
			rec_last_addr:  0,
			// "rsize": 0,
		},
		{
			// 1
			rec_start_addr: 8_065,
			seg_start_time: -120_450_514.89409208,
			seg_last_time:  510_321_600,
			int_len:        510_321_600,
			rec_last_addr:  54_071,
			// "rsize": 20,
			// init =  391_827_779.23949105
			// n =  500.0
		},
		{
			// 2
			rec_start_addr: 54_072,
			seg_start_time: -780_357_478,
			seg_last_time:  -120_450_514,
			int_len:        -120_450_514,
			rec_last_addr:  100_078,
			// "rsize": 20,
			// init =  -256344669.7390517
			// n =  500.0
		},
		{
			// 3
			rec_start_addr: 100_079,
			seg_start_time: -1_428_684_898,
			seg_last_time:  -780_357_478,
			int_len:        -780_357_478,
			rec_last_addr:  146_085,
			// "rsize": 20,
			// init =  -909271241.4267497
			// n =  500.0
		},
		{
			// 4
			rec_start_addr: 146_086,
			seg_start_time: -2_072_145_540,
			seg_last_time:  -1_428_684_898,
			int_len:        -1_428_684_898,
			rec_last_addr:  192_092,
			// "rsize": 20,
			// init =  -1559481315.7158456
			// n =  500.0
		},
		{
			// 5
			rec_start_addr: 192_093,
			seg_start_time: -2_719_741_734,
			seg_last_time:  -2_072_145_540,
			int_len:        -2_072_145_540,
			rec_last_addr:  238_099,
			// "rsize": 20,
			// init =  -2202930574.1405
			// n =  500.0
		},
		{
			// 6
			rec_start_addr: 238_100,
			seg_start_time: -3_155_716_800,
			seg_last_time:  -2_719_741_734,
			int_len:        -2_758_703_306,
			// intlen =  -2758703306.856402 ??????????????
			rec_last_addr: 268_464,
			// "rsize": 20,
			// init =  -2891244118.0133157
			// n =  330.0
		},
		{
			// 7
			rec_start_addr: 268_465,
			seg_start_time: 510_321_600,
			seg_last_time:  1_134_172_990,
			int_len:        1_134_172_990,
			rec_last_addr:  314_471,
			// "rsize": 20,
			// init =  1002901396.5747848
			// n =  500.0
		},
		{
			// 8
			rec_start_addr: 314_472,
			seg_start_time: 1_134_172_990,
			seg_last_time:  1_785_159_000,
			int_len:        1_785_159_000,
			rec_last_addr:  360_478,
			// "rsize": 20,
			// init =  1655649864.227839
			// n =  500.0
		},
		{
			// 9
			rec_start_addr: 360_479,
			seg_start_time: 1_785_159_000,
			seg_last_time:  2_451_382_230,
			int_len:        2_451_382_230,
			rec_last_addr:  406_485,
			// "rsize": 20,
			// init =  2313659513.647009
			// n =  500.0
		},
		{
			// 10
			rec_start_addr: 406_486,
			seg_start_time: 2_451_382_230,
			seg_last_time:  3_098_798_447,
			int_len:        3_098_798_447,
			rec_last_addr:  452_492,
			// "rsize": 20,
			// init =  2967630632.730112
			// n =  500.0
		},
		{
			// 11
			rec_start_addr: 452_493,
			seg_start_time: 3_098_798_447,
			seg_last_time:  3_187_252_800,
			int_len:        3_187_252_800,
			rec_last_addr:  458_842,
			// "rsize": 20,
			// init =  3185947430.156654
			// n =  69.0
		},
	}

	return calcCC(seconds, chiron, CHIRON_FILE_RECORDS)

}

func calcCeres(seconds float64, ceres []float64) types.Position {

	var CERES_FILE_RECORDS = []FileRecords{
		{
			// SUN
			rec_start_addr: 0,
			seg_start_time: 0,
			seg_last_time:  0,
			int_len:        0,
			rec_last_addr:  0,
		},
		{
			// 1
			rec_start_addr: 8_065,
			seg_start_time: 14_308_254,
			seg_last_time:  631_108_800,
			int_len:        631_108_800,
			rec_last_addr:  54_071,
		},
		{
			// 2
			rec_start_addr: 54_072,
			seg_start_time: -651_515_663,
			seg_last_time:  14_308_254,
			int_len:        14_308_254,
			rec_last_addr:  100_078,
		},
		{
			// 3
			rec_start_addr: 100_079,
			seg_start_time: -1_303_985_603,
			seg_last_time:  -651_515_663,
			int_len:        -651_515_663,
			rec_last_addr:  146_085,
		},
		{
			// 4
			rec_start_addr: 146_086,
			seg_start_time: -1_965_211_060,
			seg_last_time:  -1_303_985_603,
			int_len:        -1_303_985_603,
			rec_last_addr:  192_092,
		},
		{
			// 5
			rec_start_addr: 192_093,
			seg_start_time: -2_619_710_503,
			seg_last_time:  -1_965_211_060,
			int_len:        -1_965_211_060,
			rec_last_addr:  238_099,
		},
		{
			// 6
			rec_start_addr: 238_100,
			seg_start_time: -3_155_716_800,
			seg_last_time:  -2_619_710_503,
			int_len:        -2_638_708_105,
			rec_last_addr:  276_285,
		},
		{
			// 7
			rec_start_addr: 276_286,
			seg_start_time: -3_155_716_800,
			seg_last_time:  -2_619_710_503,
			int_len:        -2_638_708_105,
			rec_last_addr:  276_285,
			// "rsize": 20,
			// init =  1002901396.5747848
			// n =  500.0
		},
		{
			// 8
			rec_start_addr: 322_293,
			seg_start_time: 1_274_239_976,
			seg_last_time:  1_919_182_997,
			int_len:        1_919_182_997,
			rec_last_addr:  368_299,
		},
		{
			// 9
			rec_start_addr: 368_300,
			seg_start_time: 1_919_182_997,
			seg_last_time:  2_575_120_362,
			int_len:        2_575_120_362,
			rec_last_addr:  414_306,
		},
		{
			// 10
			rec_start_addr: 414_307,
			seg_start_time: 2_575_120_362,
			seg_last_time:  3_187_252_800,
			int_len:        3_092_529_717,
			rec_last_addr:  457_920,
		},
	}

	return calcCC(seconds, ceres, CERES_FILE_RECORDS)

}

func calcCC(seconds float64, pfile []float64, FILE_RECORDS []FileRecords) types.Position {

	const max_dim = 20
	const DFLSIZ = 91
	const BUFSIZ = 100
	var total_summaries_number = len(FILE_RECORDS) - 1
	// fmt.Println("total_summaries_number = ", total_summaries_number)

	for i_summ := 1; i_summ <= total_summaries_number; i_summ++ {
		if FILE_RECORDS[i_summ].seg_start_time < seconds &&
			FILE_RECORDS[i_summ].seg_last_time > seconds {

			start_adress := FILE_RECORDS[i_summ].rec_start_addr
			last_adress := FILE_RECORDS[i_summ].rec_last_addr

			n_of_rec := int(pfile[last_adress-1])

			// Number of directory epochs
			var n_of_dir int = n_of_rec / BUFSIZ

			OFFD := last_adress - n_of_dir - 2
			OFFE := OFFD - n_of_rec

			RECNO := -1

			if n_of_rec <= BUFSIZ {

				data := pfile[(OFFE):(OFFE + n_of_rec)]

				if seconds < data[0] || seconds > data[len(data)-1] {
					fmt.Println("we have a problem")
				}

				for i, v := range data {
					if v == seconds {
						// fmt.Println("equality , index =", i, "value = ", arr[i])
						RECNO = i
						break
					}
					if v > seconds {
						// fmt.Println("index =", i-1, "value = ", arr[i-1])
						RECNO = i - 1
						break
					}
				}

			} else {
				for dir := 0; dir < n_of_dir; dir++ {

					data := pfile[(OFFD + dir) : (OFFD+dir)+1][0]

					if data > seconds {
						OFFD = OFFE + (dir)*BUFSIZ
						data := pfile[(OFFD):(OFFD + BUFSIZ)]

						// looking for the largest array element less than  dateInSeconds
						// and get its index
						if seconds < data[0] || seconds > data[len(data)-1] {
							fmt.Println("we have a problem")
						}

						for i, v := range data {
							if v == seconds {
								// fmt.Println("equality , index =", i, "value = ", arr[i])
								RECNO = i
								break
							}
							if v > seconds {
								// fmt.Println("index =", i-1, "value = ", arr[i-1])
								RECNO = i - 1
								break
							}
						}

						RECNO = dir*BUFSIZ + RECNO
						break
					}
				}
			}

			if RECNO == -1 {
				// print("chiron final records sec = ", dateInSeconds)
				Ind := n_of_rec % BUFSIZ
				data := pfile[(last_adress - n_of_dir - Ind) : last_adress-n_of_dir]

				// looking for the largest array element less than  dateInSeconds
				// and get its index
				if seconds < data[0] || seconds > data[len(data)-1] {
					fmt.Println("we have a problem")
				}

				for i, v := range data {
					if v == seconds {
						// fmt.Println("equality , index =", i, "value = ", arr[i])
						RECNO = i
						break
					}
					if v > seconds {
						// fmt.Println("index =", i-1, "value = ", arr[i-1])
						RECNO = i - 1
						break
					}
				}

				RECNO = (n_of_dir)*BUFSIZ + RECNO
			}

			OFFR := start_adress - 1 + (RECNO)*DFLSIZ
			mda_record := pfile[(OFFR):(OFFR + DFLSIZ)]

			// print("dateInSeconds = ", dateInSeconds)
			// print(mda_record)

			TL := mda_record[0]
			G := mda_record[1 : max_dim+1]
			REFPOS := []float64{mda_record[max_dim+1], mda_record[max_dim+3], mda_record[max_dim+5]}
			REFVEL := []float64{mda_record[max_dim+2], mda_record[max_dim+4], mda_record[max_dim+6]}

			KQMAX1 := int(mda_record[4*max_dim+7])

			KQ := []float64{(mda_record[4*max_dim+8]), (mda_record[4*max_dim+9]), (mda_record[4*max_dim+10])}

			KS := KQMAX1 - 1
			MQ2 := KQMAX1 - 2
			DELTA := seconds - TL
			TP := DELTA

			// FC = [0 for i in range(max_dim)]
			var FC [max_dim]float64
			FC[0] = 1.0
			// WC = [0 for i in range(max_dim - 1)]
			var WC [max_dim - 1]float64

			for J := 1; J < MQ2+1; J++ {
				if G[J-1] == 0.0 {
					fmt.Println("SPKE21\nA value of zero was found at index {0} of the step size vector.")
				}
				FC[J] = TP / G[J-1]
				WC[J-1] = DELTA / G[J-1]
				TP = DELTA + G[J-1]
			}
			// W = [0 for i in range(max_dim + 2)]
			var W [max_dim + 2]float64

			//
			//     Collect KQMAX1 reciprocals.
			//     KS = KQMAX1 - 1     KS = KQMAX1 - 1
			// for J in range(1, KQMAX1 + 1):
			//     W[J - 1] = 1.0 / float(J)
			for J := 1; J < KQMAX1+1; J++ {
				W[J-1] = 1.0 / float64(J)
			}

			//
			//     Compute the W(K) terms needed for the position interpolation
			//     (Note,  it is assumed throughout this routine that KS, which
			//     starts out as KQMAX1-1 (the ``maximum integration'')
			//     is at least 2.
			//

			JX := 0
			KS1 := KS - 1

			for KS >= 2 {
				JX = JX + 1

				for J := 1; J < JX+1; J++ {
					W[J+KS-1] = FC[J]*W[J+KS1-1] - WC[J-1]*W[J+KS-1]
				}
				KS = KS1
				KS1 = KS1 - 1
			}
			//
			//     Perform position interpolation: (Note that KS = 1 right now.
			//     We don't know much more than that.)
			//

			STATE := []float64{0, 0, 0}

			// DTtest = np.reshape(  mda_record[max_dim + 7 : max_dim * 4 + 7], (max_dim, 3), order="F"  )
			first_arr := mda_record[max_dim+7 : max_dim*2+7]
			second_arr := mda_record[max_dim*2+7 : max_dim*3+7]
			third_arr := mda_record[max_dim*3+7 : max_dim*4+7]
			DTtest := [][]float64{first_arr, second_arr, third_arr}

			for Ii := 0; Ii < 3; Ii++ {
				KQQ := KQ[Ii]
				SUM := 0.0

				for J := KQQ; J > 0; J-- {
					// v SUM = SUM + DTtest[J - 1][Ii] * W[J - 1 + KS]
					// SUM = SUM + DTtest[int(J-1)][Ii]*W[(int(J)-1+KS)]
					SUM = SUM + DTtest[Ii][int(J-1)]*W[(int(J)-1+KS)]
				}

				STATE[Ii] = REFPOS[Ii] + DELTA*(REFVEL[Ii]+DELTA*SUM)
			}
			return types.Position{X: STATE[0], Y: STATE[1], Z: STATE[2]}
		}

	}

	return types.Position{X: 0, Y: 0, Z: 0}
}
