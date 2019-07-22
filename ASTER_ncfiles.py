#go into the main location directory
#list=`ls */*.nc` ; for item in $list ; do echo \'$item\' \,; done

dungeness_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/dungeness'
dungeness_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/dungeness_1999/nighttime'

hartlepool_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hartlepool'
hartlepool_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hartlepool_1999/nighttime'

heysham_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/heysham'
heysham_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/heysham_1999/nighttime'

hinkley_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hinkley'
hinkley_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hinkley_1999/nighttime'

hunterston_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hunterston/nighttime'
hunterston_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/hunterston_1999/nighttime'

sizewell_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/sizewell'
sizewell_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/sizewell_1999/nighttime'

torness_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/torness'
torness_1999_path = '/storage/silver/surft/users/mp877190/data/datastore/EE/ASTER_L1T/torness_1999/nighttime'

dungeness_tir = ['2176209107_tir/AST_L1T_00307142015213909_20150715182923_26770.nc' ,
'2201043071_tir/AST_L1T_00311012015215052_20151104012453_28462.nc' ,
'2201043137_tir/AST_L1T_00311012015215043_20151104012443_28414.nc' ,
'2204764632_tir/AST_L1T_00312122015214426_20151215191616_10742.nc' ,
'2244625735_tir/AST_L1T_00312052016215021_20161206085615_17459.nc' ,
'2244625738_tir/AST_L1T_00312052016215030_20161206085625_18673.nc' ,
'2244682320_tir/AST_L1T_00312072016213817_20161208153635_5890.nc' ,
'2244682321_tir/AST_L1T_00312072016213808_20161208153545_5038.nc' ,
'2246988030_tir/AST_L1T_00302232017215004_20170224111543_25411.nc' ,
'2246988070_tir/AST_L1T_00302232017215013_20170224111553_25437.nc' ,
'2248725898_tir/AST_L1T_00304122017215007_20170413081452_2614.nc' ,
'2250350727_tir/AST_L1T_00305072017214412_20170508094440_4500.nc' ,
#'2251707868_tir/AST_L1T_00305142017215023_20170515112318_23809.nc' ,
'2251707875_tir/AST_L1T_00305142017215014_20170515112308_23605.nc' ,
'2273056952_tir/AST_L1T_00308182017215045_20170819112220_8676.nc' ,
'2273057100_tir/AST_L1T_00308182017215036_20170819114217_15136.nc' ,
'2273088987_tir/AST_L1T_00308202017213824_20170821105036_11005.nc' ,
'2273088988_tir/AST_L1T_00308202017213832_20170821105046_12130.nc']

dungeness_all = ['2190893700/AST_L1T_00309102015110443_20150917155627_11321.nc' ,
'2190908449/AST_L1T_00309102015110434_20150917170759_9563.nc']

dungeness_1999_tir = ['2176209107/AST_L1T_00307142015213909_20150715182923_2677.nc' ,
'2201043071-useless/AST_L1T_00311012015215052_20151104012453_2846.nc' ,
'2201043137-maybe/AST_L1T_00311012015215043_20151104012443_2841.nc' ,
'2204764632-useless/AST_L1T_00312122015214426_20151215191616_1074.nc' ,
'2244625735/AST_L1T_00312052016215021_20161206085615_1745.nc' ,
'2244625738-useless/AST_L1T_00312052016215030_20161206085625_1867.nc' ,
'2244682321/AST_L1T_00312072016213808_20161208153545_5038.nc' ,
'2246988030-maybe/AST_L1T_00302232017215004_20170224111543_2541.nc' ,
'2246988070-useless/AST_L1T_00302232017215013_20170224111553_2543.nc' ,
'2248725898-useless/AST_L1T_00304122017215007_20170413081452_2614.nc' ,
'2250350727-maybe/AST_L1T_00305072017214412_20170508094440_4500.nc' ,
#'2251707868-useless/AST_L1T_00305142017215023_20170515112318_2380.nc' ,
'2251707875-maybe/AST_L1T_00305142017215014_20170515112308_2360.nc' ,
'2273056952-useless/AST_L1T_00308182017215045_20170819112220_8676.nc' ,
'2273057100-maybe/AST_L1T_00308182017215036_20170819114217_1513.nc' ,
'2273088987/AST_L1T_00308202017213824_20170821105036_1100.nc' ,
'2282050530/AST_L1T_00311242017213807_20171128015051_2473.nc' ,
'2283338346-useless/AST_L1T_00312172017214410_20171218101922_2075.nc' ,
'2283848645-useless/AST_L1T_00312242017215023_20171225094521_1113.nc' ,
'2283848654/AST_L1T_00312242017215014_20171225094711_1764.nc']

hartlepool_tir = ['2170752441_tir/AST_L1T_2170752441.nc' ,
'2176209612_tir/AST_L1T_00307142015214019_20150715183143_28559.nc' ,
#'2176209614_tir/AST_L1T_00307142015214011_20150715183133_28553.nc' ,
'2179288444_tir/AST_L1T_00307282015215232_20150729140650_6809.nc' ,
'2179288449_tir/AST_L1T_00307282015215223_20150729140520_2273.nc' ,
'2181374236_tir/AST_L1T_00308062015214619_20150807183925_24828.nc' ,
'2204764631_tir/AST_L1T_00312122015214528_20151215191806_12319.nc' ,
'2212902755_tir/AST_L1T_00302122016215744_20160216034903_4135.nc' ,
'2242834653_tir/AST_L1T_00310022016215133_20161003112949_25905.nc' ,
'2242834655_tir/AST_L1T_00310022016215124_20161003112909_25142.nc' ,
'2244625737_tir/AST_L1T_00312052016215123_20161206085525_14901.nc' ,
#'2244625739_tir/AST_L1T_00312052016215132_20161206085535_15093.nc' ,
#'2244682329_tir/AST_L1T_00312072016213919_20161208153835_7664.nc' ,
#'2244682331_tir/AST_L1T_00312072016213910_20161208153755_6824.nc' ,
#'2246194802_tir/AST_L1T_00301292017215718_20170130093752_6854.nc' ,
#'2246203709_tir/AST_L1T_00301292017215727_20170130124231_5023.nc' ,
#'2250371768_tir/AST_L1T_00305072017214514_20170508121521_30028.nc' ,
#'2263075900_tir/AST_L1T_00306242017214515_20170625165933_17314.nc' ,
'2264743330_tir/AST_L1T_00307082017215739_20170709094314_4113.nc' ,
#'2264743422_tir/AST_L1T_00307082017215730_20170709094815_19954.nc' ,
'2267510657_tir/AST_L1T_00307242017215739_20170725145250_328.nc' ,
'2269788158_tir/AST_L1T_00307012017215116_20170804225828_501.nc' ,
'2273056946_tir/AST_L1T_00308182017215138_20170819112340_12688.nc' ,
'2273057071_tir/AST_L1T_00308182017215147_20170819113727_32223.nc' ,
'2273088872_tir/AST_L1T_00308202017213934_20170821105006_9030.nc' ,
#'2273088874_tir/AST_L1T_00308202017213925_20170821104926_7705.nc' ,
'2273148673_tir/AST_L1T_00308092017215746_20170822150713_23300.nc']

hartlepool_all = [#'2169414789/2169414789.nc' ,
'2189214848/AST_L1T_00309062015112811_20150908212613_18107.nc' ,
'2189215047/AST_L1T_00309062015112820_20150908212743_20315.nc' ,
'2266466308/AST_L1T_00307182017112138_20170719134949_22714.nc']

hartlepool_1999_tir = ['2167794869-maybe/AST_L1T_00312172011214458_20150608183653_5534.nc' ,
'2170752441/AST_L1T_00303122014214502_20150620065510_7081.nc' ,
'2176209612-maybe/AST_L1T_00307142015214019_20150715183143_2855.nc' ,
'2176209614-useless/AST_L1T_00307142015214011_20150715183133_2855.nc' ,
'2179288444-maybe/AST_L1T_00307282015215232_20150729140650_6809.nc' ,
'2181374236-useless/AST_L1T_00308062015214619_20150807183925_2482.nc' ,
'2204764631/AST_L1T_00312122015214528_20151215191806_1231.nc' ,
'2212902755-useless/AST_L1T_00302122016215744_20160216034903_4135.nc' ,
'2242834653/AST_L1T_00310022016215133_20161003112949_2590.nc' ,
'2244625737-useless/AST_L1T_00312052016215123_20161206085525_1490.nc' ,
'2244625739/AST_L1T_00312052016215132_20161206085535_1509.nc' ,
'2244682329-useless/AST_L1T_00312072016213919_20161208153835_7664.nc' ,
'2246194802-useless/AST_L1T_00301292017215718_20170130093752_6854.nc' ,
'2246203709-useless/AST_L1T_00301292017215727_20170130124231_5023.nc' ,
'2250371768/AST_L1T_00305072017214514_20170508121521_3002.nc' ,
'2263075900/AST_L1T_00306242017214515_20170625165933_1731.nc' ,
'2264743422-useless/AST_L1T_00307082017215730_20170709094815_1995.nc' ,
'2267510657-useless/AST_L1T_00307242017215739_20170725145250_328..nc' ,
'2273056946-useless/AST_L1T_00308182017215138_20170819112340_1268.nc' ,
'2273057071/AST_L1T_00308182017215147_20170819113727_3222.nc' ,
'2273088872-maybe/AST_L1T_00308202017213934_20170821105006_9030.nc' ,
#'2273088874-useless/AST_L1T_00308202017213925_20170821104926_7705.nc' ,
'2273148673-useless/AST_L1T_00308092017215746_20170822150713_2330.nc' ,
'2283268940-useless/AST_L1T_00312152017215733_20171216104427_7116.nc' ,
'2283268965-useless/AST_L1T_00312152017215724_20171216104347_6204.nc' ,
'2283327817/AST_L1T_00312172017214512_20171218090011_1072.nc']

heysham_tir = [#ebb,
 '2194252766_tir/AST_L1T_00309302015215121_20151001103139_9227.nc' ] #,
"""
#flood,
 '2195766298_tir/AST_L1T_00310072015215738_20151008173904_4582.nc' ,
'2195766302_tir/AST_L1T_00310072015215729_20151008173854_4565.nc' ,
'2204905728_tir/AST_L1T_00312102015215743_20151216233925_22227.nc' ,
'2204905729_tir/AST_L1T_00312102015215734_20151216233915_22187.nc' ,
'2215166184_tir/AST_L1T_00303062016220321_20160307094533_30663.nc' ,
#flood,
 '2225511855_tir/AST_L1T_00305092016220357_20160510092009_24656.nc' ,
'2237715318_tir/AST_L1T_00307302016215200_20160731223902_3071.nc' ,
'2243406412_tir/AST_L1T_00310252016215722_20161026141904_29965.nc' ,
'2243411732_tir/AST_L1T_00310252016215731_20161026170436_3631.nc' ,
'2244972831_tir/AST_L1T_00312192016220334_20161220130946_21952.nc' ,
#flood,
 '2280337814_tir/AST_L1T_00310052017215132_20171006104449_4958.nc']
"""
heysham_all = ['2171270613/AST_L1T_00307242014113352_20150622072320_92683.nc']

heysham_1999_tir = [#ebb'
 '2194252766/AST_L1T_00309302015215121_20151001103139_9227.nc' ,
#flood,
# '2195766298/AST_L1T_00310072015215738_20151008173904_4582.nc' ,
#'2204905728-maybe/AST_L1T_00312102015215743_20151216233925_2222.nc' ,
'2215166184-maybe/AST_L1T_00303062016220321_20160307094533_3066.nc' ,
#'2225511855-maybe/AST_L1T_00305092016220357_20160510092009_2465.nc' ,
#ebb,
 '2237715318-useless/AST_L1T_00307302016215200_20160731223902_3071.nc' ,
'2243411732/AST_L1T_00310252016215731_20161026170436_3631.nc' ,
'2244972831/AST_L1T_00312192016220334_20161220130946_2195.nc' ]
#flood,
# '2280337814/AST_L1T_00310052017215132_20171006104449_4958.nc']

hinkley_tir = []
hinkley_all = []

hinkley_1999_tir = ['2196532773/AST_L1T_00310142015220302_20151015094600_3066.nc' ,
'2204431959/AST_L1T_00312082015220914_20151209160827_2693.nc' ,
'2227419112/AST_L1T_00305162016220930_20160517112825_319..nc' ,
'2227974100-useless/AST_L1T_00305182016215719_20160519234652_2438.nc' ,
'2244544878/AST_L1T_00312032016220251_20161204124221_1904.nc' ,
'2246717037-useless/AST_L1T_00302142017215623_20170215131445_3219.nc' ,
'2246717053-maybe/AST_L1T_00302142017215632_20170215131455_3251.nc' ,
'2276747697-useless/AST_L1T_00309262017215703_20170927133503_2525.nc' ,
'2281828159-maybe/AST_L1T_00311132017215650_20171114091503_1712.nc' ,
'2284203502-maybe/AST_L1T_00301142018220857_20180115094448_2935.nc']

hunterston_tir = ['2194253176_tir_useless/AST_L1T_00309302015215156_20151001103249_1016.nc' ,
'2196533541_tir/AST_L1T_00310142015220422_20151015094850_911..nc' ,
'2205486063_tir/AST_L1T_00312262015215810_20151229092219_7115.nc' ,
'2206378064_tir/AST_L1T_00301182016220429_20160119081029_7195.nc' ,
'2243173336_tir/AST_L1T_00310162016220410_20161017113603_1089.nc' ,
'2243668502_tir/AST_L1T_00311012016220408_20161102132849_9262.nc' ,
'2246717299_tir/AST_L1T_00302142017215743_20170215072224_2277.nc' ,
'2280336377_tir/AST_L1T_00310052017215159_20171006104319_3232.nc' ,
'2280337818_tir/AST_L1T_00310052017215208_20171006104519_6154.nc' ,
'2281306019_tir_useless/AST_L1T_00310282017215800_20171029105550_2396.nc']

hunterston_all = []

hunterston_1999_tir = ['2194253176/AST_L1T_00309302015215156_20151001103249_1016.nc' ,
'2196533541/AST_L1T_00310142015220422_20151015094850_911..nc' ,
'2205486063/AST_L1T_00312262015215810_20151229092219_7115.nc' ,
'2206378064/AST_L1T_00301182016220429_20160119081029_7195.nc' ,
'2243173336/AST_L1T_00310162016220410_20161017113603_1089.nc' ,
'2243668502/AST_L1T_00311012016220408_20161102132849_9262.nc' ,
'2246717299/AST_L1T_00302142017215743_20170215072224_2277.nc' ,
'2280336377/AST_L1T_00310052017215159_20171006104319_3232.nc' ,
'2280337818/AST_L1T_00310052017215208_20171006104519_6154.nc' ,
'2281306019/AST_L1T_00310282017215800_20171029105550_2396.nc' ,
'2282091413/AST_L1T_00311292017215757_20171130100953_4747.nc' ,
'2284203493/AST_L1T_00301142018221008_20180115094108_1687.nc' ,
'2284203495/AST_L1T_00301142018221017_20180115094118_1713.nc']

sizewell_tir = ['2171658987_tir/AST_L1T_00306192015214523_20150623235018_16073.nc' ,
'2212460343_tir/AST_L1T_00302072016213848_20160212042549_16748.nc' ,
'2212460344_tir/AST_L1T_00302072016213839_20160212042519_14888.nc' ,
'2243084601_tir/AST_L1T_00310112016214437_20161012175058_7260.nc' ,
'2243301057_tir/AST_L1T_00310202016213823_20161021083319_2349.nc' ,
'2243303277_tir/AST_L1T_00310202016213832_20161021105651_21137.nc' ,
'2245528172_tir/AST_L1T_00301082017213825_20170110145958_12795.nc' ,
'2245528176_tir/AST_L1T_00301082017213816_20170110145948_12770.nc' ,
'2248406074_tir/AST_L1T_00303292017213814_20170330094603_20980.nc' ,
'2248406076_tir/AST_L1T_00303292017213823_20170330094613_21650.nc' ,
'2254262702_tir/AST_L1T_00305252017213216_20170526123109_6223.nc' ,
'2281127042_tir/AST_L1T_00310232017213826_20171024105853_12245.nc' ,
'2281127053_tir/AST_L1T_00310232017213835_20171024110003_15513.nc']

sizewell_all = ['2169538372/AST_L1T_00306072013110958_20150616003108_95527.nc' ,
'2169538377/AST_L1T_00306072013110949_20150616003108_95529.nc' ,
'2190909118/AST_L1T_00309102015110416_20150917171151_23916.nc']

sizewell_1999_tir = ['2161044357-useless/AST_L1T_00302142007213837_20150518070934_8937.nc' ,
'2171658987-maybe/AST_L1T_00306192015214523_20150623235018_1607.nc' ,
'2212460343-useless/AST_L1T_00302072016213848_20160212042549_1674.nc' ,
'2212460344-useless/AST_L1T_00302072016213839_20160212042519_1488.nc' ,
'2243084601-useless/AST_L1T_00310112016214437_20161012175058_7260.nc' ,
'2243301057-maybe/AST_L1T_00310202016213823_20161021083319_2349.nc' ,
'2243303277-useless/AST_L1T_00310202016213832_20161021105651_2113.nc' ,
'2245528172-useless/AST_L1T_00301082017213825_20170110145958_1279.nc' ,
'2245528176-useless/AST_L1T_00301082017213816_20170110145948_1277.nc' ,
'2248406074-maybe/AST_L1T_00303292017213814_20170330094603_2098.nc' ,
'2248406076-useless/AST_L1T_00303292017213823_20170330094613_2165.nc' ,
'2254262702-useless/AST_L1T_00305252017213216_20170526123109_6223.nc' ,
'2281127042-maybe/AST_L1T_00310232017213826_20171024105853_1224.nc' ,
'2281127053-useless/AST_L1T_00310232017213835_20171024110003_1551.nc' ,
'2284162482-maybe/AST_L1T_00301112018213821_20180112093242_4784.nc' ,
'2284162561-useless/AST_L1T_00301112018213830_20180112094933_4742.nc' ,
'2284385061-useless/AST_L1T_00301182018214431_20180119101454_1166.nc']

torness_tir = ['2179312100_tir/AST_L1T_00307282015215258_20150729164916_15014.nc' ,
'2179312102_tir/AST_L1T_00307282015215250_20150729164906_13786.nc' ,
'2181342046_tir/AST_L1T_00308062015214645_20150807153023_532.nc' ,
'2203578426_tir/AST_L1T_00311242015215822_20151125224310_27884.nc' ,
'2212902849_tir/AST_L1T_00302122016215811_20160216034953_5045.nc' ,
'2223904640_tir/AST_L1T_00305022016215815_20160503112108_22698.nc' ,
'2242834656_tir/AST_L1T_00310022016215151_20161003113009_26102.nc' ,
'2242834730_tir/AST_L1T_00310022016215200_20161003113049_26853.nc' ,
'2244267224_tir/AST_L1T_00311282016214547_20161129083103_30595.nc' ,
'2244777569_tir/AST_L1T_00312122016215800_20161213090933_26016.nc' ,
'2246203708_tir/AST_L1T_00301292017215745_20170130124251_5184.nc' ,
'2263071830_tir/AST_L1T_00306242017214541_20170625135407_20458.nc' ,
'2264743284_tir/AST_L1T_00307082017215757_20170709094051_28423.nc' ,
'2267510594_tir/AST_L1T_00307242017215806_20170725144839_26204.nc' ,
'2269788156_tir/AST_L1T_00307012017215143_20170804225717_30049.nc']
#'2273148667_tir/AST_L1T_00308092017215813_20170822150803_24367.nc']

torness_all = []

torness_1999_tir = ['2179312100-useless/AST_L1T_00307282015215258_20150729164916_1501.nc' ,
'2179312102-maybe/AST_L1T_00307282015215250_20150729164906_1378.nc' ,
'2181342046/AST_L1T_00308062015214645_20150807153023_532..nc' ,
'2203578426-useless/AST_L1T_00311242015215822_20151125224310_2788.nc' ,
'2212902849-maybe/AST_L1T_00302122016215811_20160216034953_5045.nc' ,
'2223904640-useless/AST_L1T_00305022016215815_20160503112108_2269.nc' ,
'2242834656/AST_L1T_00310022016215151_20161003113009_2610.nc' ,
'2244267224/AST_L1T_00311282016214547_20161129083103_3059.nc' ,
'2244777569-useless/AST_L1T_00312122016215800_20161213090933_2601.nc' ,
'2246203708/AST_L1T_00301292017215745_20170130124251_5184.nc' ,
'2263071830-maybe/AST_L1T_00306242017214541_20170625135407_2045.nc' ,
'2264743284/AST_L1T_00307082017215757_20170709094051_2842.nc' ,
'2267510594-maybe/AST_L1T_00307242017215806_20170725144839_2620.nc' ,
'2269788156/AST_L1T_00307012017215143_20170804225717_3004.nc' ,
'2273148667/AST_L1T_00308092017215813_20170822150803_2436.nc' ,
'2283268487/AST_L1T_00312152017215751_20171216101135_1491.nc' ,
'2283327816/AST_L1T_00312172017214539_20171218090131_1591.nc']
