#!/usr/bin/python

small_s = ["EG_transcript_15350", "EG_transcript_12671", "EG_transcript_19007", "EG_transcript_14195", "EG_transcript_16866", "EG_transcript_19007", "EG_transcript_14195", "EG_transcript_16866", "EG_transcript_19007", "EG_transcript_14195", "EG_transcript_16866", "EG_transcript_54525", "EG_transcript_18154", "EG_transcript_53084", "EG_transcript_56241", "EG_transcript_35367", "EG_transcript_52542", "EG_transcript_3176", "EG_transcript_24541", "EG_transcript_13997", "EG_transcript_3042", "EG_transcript_7178", "EG_transcript_7444", "EG_transcript_7178", "EG_transcript_7444", "EG_transcript_43242", "EG_transcript_5584", "EG_transcript_23105", "EG_transcript_790", "EG_transcript_18811", "EG_transcript_26745", "EG_transcript_443", "EG_transcript_7444", "EG_transcript_7178", "EG_transcript_26745", "EG_transcript_21240", "EG_transcript_23105", "EG_transcript_33577", "EG_transcript_15718", "EG_transcript_23105", "EG_transcript_12637", "EG_transcript_39982", "EG_transcript_44308", "EG_transcript_2112", "EG_transcript_15810", "EG_transcript_21954", "EG_transcript_30136", "EG_transcript_21954", "EG_transcript_30136", "EG_transcript_20817", "EG_transcript_20817", "EG_transcript_32001", "EG_transcript_40131", "EG_transcript_11597", "EG_transcript_21706", "EG_transcript_25501", "EG_transcript_40456", "EG_transcript_23105", "EG_transcript_17968", "EG_transcript_52320", "EG_transcript_790", "EG_transcript_4714", "EG_transcript_4714", "EG_transcript_21511", "EG_transcript_26745", "EG_transcript_23963", "EG_transcript_21240", "EG_transcript_23963", "EG_transcript_24888", "EG_transcript_36591", "EG_transcript_30657", "EG_transcript_17529", "EG_transcript_12671", "EG_transcript_6467", "EG_transcript_26745", "EG_transcript_21240", "EG_transcript_17106", "EG_transcript_15718", "EG_transcript_23844", "EG_transcript_24567", "EG_transcript_25719"]
big_s = ["EG_transcript_9781", "EG_transcript_8016", "EG_transcript_2280", "EG_transcript_28774", "EG_transcript_8907", "EG_transcript_3799", "EG_transcript_20319", "EG_transcript_16283", "EG_transcript_8914", "EG_transcript_21219", "EG_transcript_1993", "EG_transcript_6291", "EG_transcript_2696", "EG_transcript_3547", "EG_transcript_14917", "EG_transcript_1177", "EG_transcript_11341", "EG_transcript_33313", "EG_transcript_12260", "EG_transcript_1196", "EG_transcript_9064", "EG_transcript_4750", "EG_transcript_12306", "EG_transcript_10940", "EG_transcript_1290", "EG_transcript_1536", "EG_transcript_12480", "EG_transcript_1290", "EG_transcript_1245", "EG_transcript_5", "EG_transcript_8289", "EG_transcript_9695", "EG_transcript_15919", "EG_transcript_5917", "EG_transcript_97", "EG_transcript_44308", "EG_transcript_52992", "EG_transcript_30159", "EG_transcript_8679", "EG_transcript_19329", "EG_transcript_16724", "EG_transcript_13352", "EG_transcript_3774", "EG_transcript_4713", "EG_transcript_1241", "EG_transcript_8320", "EG_transcript_1245", "EG_transcript_19235", "EG_transcript_18332", "EG_transcript_8471", "EG_transcript_20567", "EG_transcript_18179", "EG_transcript_1548", "EG_transcript_7356", "EG_transcript_11577", "EG_transcript_2549", "EG_transcript_7357", "EG_transcript_7357", "EG_transcript_7357", "EG_transcript_7357", "EG_transcript_7357", "EG_transcript_7967", "EG_transcript_40410", "EG_transcript_3042", "EG_transcript_2608", "EG_transcript_1951", "EG_transcript_7444", "EG_transcript_7444", "EG_transcript_1241", "EG_transcript_8047", "EG_transcript_7826", "EG_transcript_27641", "EG_transcript_1163", "EG_transcript_2741", "EG_transcript_4771", "EG_transcript_7", "EG_transcript_4764", "EG_transcript_4905", "EG_transcript_11779", "EG_transcript_975", "EG_transcript_19995", "EG_transcript_10192", "EG_transcript_12357", "EG_transcript_4934", "EG_transcript_21796", "EG_transcript_908", "EG_transcript_2494", "EG_transcript_6076", "EG_transcript_15319", "EG_transcript_5258", "EG_transcript_27059", "EG_transcript_24213", "EG_transcript_29716", "EG_transcript_12332", "EG_transcript_11456", "EG_transcript_9664", "EG_transcript_7994", "EG_transcript_4731", "EG_transcript_3639", "EG_transcript_17220", "EG_transcript_16406", "EG_transcript_12795", "EG_transcript_5917", "EG_transcript_5917", "EG_transcript_5917", "EG_transcript_3075", "EG_transcript_61646", "EG_transcript_668", "EG_transcript_9682", "EG_transcript_1653", "EG_transcript_27366", "EG_transcript_181", "EG_transcript_25242", "EG_transcript_6941", "EG_transcript_14945", "EG_transcript_34325", "EG_transcript_14805", "EG_transcript_5366", "EG_transcript_8362", "EG_transcript_10635", "EG_transcript_20956", "EG_transcript_14931", "EG_transcript_8425", "EG_transcript_8425", "EG_transcript_791", "EG_transcript_3371", "EG_transcript_3371", "EG_transcript_3371", "EG_transcript_12970", "EG_transcript_12621", "EG_transcript_21797", "EG_transcript_24958", "EG_transcript_3882", "EG_transcript_403", "EG_transcript_6017", "EG_transcript_5257", "EG_transcript_17132", "EG_transcript_17535", "EG_transcript_6968", "EG_transcript_9899", "EG_transcript_10083", "EG_transcript_21118", "EG_transcript_5575", "EG_transcript_5750", "EG_transcript_7178", "EG_transcript_7178", "EG_transcript_1783", "EG_transcript_12376", "EG_transcript_8986", "EG_transcript_21720", "EG_transcript_23666", "EG_transcript_14010", "EG_transcript_10850", "EG_transcript_26871", "EG_transcript_8975", "EG_transcript_924", "EG_transcript_17460", "EG_transcript_18594", "EG_transcript_2666", "EG_transcript_15154", "EG_transcript_5458", "EG_transcript_5458", "EG_transcript_5458", "EG_transcript_1376", "EG_transcript_9827", "EG_transcript_19961", "EG_transcript_7826", "EG_transcript_11927", "EG_transcript_6030", "EG_transcript_5917", "EG_transcript_10100", "EG_transcript_5076", "EG_transcript_8779", "EG_transcript_13724", "EG_transcript_10393", "EG_transcript_6017", "EG_transcript_14193", "EG_transcript_6467", "EG_transcript_5458", "EG_transcript_3017", "EG_transcript_43548", "EG_transcript_966", "EG_transcript_11592", "EG_transcript_10580", "EG_transcript_18631", "EG_transcript_11597", "EG_transcript_21143", "EG_transcript_29186", "EG_transcript_1330", "EG_transcript_27131", "EG_transcript_9064", "EG_transcript_19566", "EG_transcript_6270", "EG_transcript_3655", "EG_transcript_3760", "EG_transcript_2884", "EG_transcript_14837", "EG_transcript_3760", "EG_transcript_23774", "EG_transcript_403", "EG_transcript_3655", "EG_transcript_6933", "EG_transcript_4714", "EG_transcript_4714", "EG_transcript_12086", "EG_transcript_8517", "EG_transcript_8722", "EG_transcript_7940", "EG_transcript_11845", "EG_transcript_2393", "EG_transcript_21444", "EG_transcript_15855", "EG_transcript_3342", "EG_transcript_40456", "EG_transcript_18452", "EG_transcript_8385", "EG_transcript_36591", "EG_transcript_5437", "EG_transcript_9159", "EG_transcript_9159", "EG_transcript_11118", "EG_transcript_14424", "EG_transcript_3060", "EG_transcript_3655", "EG_transcript_3655", "EG_transcript_3655", "EG_transcript_9338", "EG_transcript_19007", "EG_transcript_19007", "EG_transcript_19007", "EG_transcript_6313", "EG_transcript_5157", "EG_transcript_5157", "EG_transcript_5157", "EG_transcript_3387", "EG_transcript_11568", "EG_transcript_10824", "EG_transcript_21271", "EG_transcript_3632", "EG_transcript_3655", "EG_transcript_8223", "EG_transcript_10131", "EG_transcript_15185", "EG_transcript_31803", "EG_transcript_10397", "EG_transcript_6432", "EG_transcript_21511", "EG_transcript_9355", "EG_transcript_4275", "EG_transcript_9425", "EG_transcript_15209", "EG_transcript_25308", "EG_transcript_10283", "EG_transcript_15475", "EG_transcript_8132", "EG_transcript_646", "EG_transcript_16382", "EG_transcript_646", "EG_transcript_4466", "EG_transcript_451", "EG_transcript_15674", "EG_transcript_485", "EG_transcript_41111", "EG_transcript_14139", "EG_transcript_485", "EG_transcript_17164", "EG_transcript_7178", "EG_transcript_21350", "EG_transcript_9632", "EG_transcript_26116", "EG_transcript_26116", "EG_transcript_16558", "EG_transcript_4398", "EG_transcript_16185", "EG_transcript_16185", "EG_transcript_12862", "EG_transcript_2633", "EG_transcript_2633", "EG_transcript_6368", "EG_transcript_7444", "EG_transcript_10250", "EG_transcript_14055", "EG_transcript_24405", "EG_transcript_18811", "EG_transcript_4109", "EG_transcript_8021", "EG_transcript_10397", "EG_transcript_2565", "EG_transcript_13767", "EG_transcript_3806", "EG_transcript_16896", "EG_transcript_30512", "EG_transcript_5361", "EG_transcript_13224", "EG_transcript_7532", "EG_transcript_14495", "EG_transcript_950", "EG_transcript_1121", "EG_transcript_9640", "EG_transcript_10182", "EG_transcript_7", "EG_transcript_10476", "EG_transcript_8825", "EG_transcript_310", "EG_transcript_14896", "EG_transcript_12587", "EG_transcript_14387", "EG_transcript_5925", "EG_transcript_1707", "EG_transcript_8197", "EG_transcript_14769", "EG_transcript_38453", "EG_transcript_15592", "EG_transcript_5", "EG_transcript_34143", "EG_transcript_15240", "EG_transcript_15887", "EG_transcript_1083", "EG_transcript_8430", "EG_transcript_15039", "EG_transcript_2037", "EG_transcript_5944", "EG_transcript_451", "EG_transcript_10864", "EG_transcript_7879", "EG_transcript_13548", "EG_transcript_19861", "EG_transcript_4447", "EG_transcript_20340", "EG_transcript_5683", "EG_transcript_276", "EG_transcript_25374", "EG_transcript_14503", "EG_transcript_21133", "EG_transcript_26209", "EG_transcript_1542", "EG_transcript_1060", "EG_transcript_4778", "EG_transcript_19736", "EG_transcript_8430", "EG_transcript_2040", "EG_transcript_1497", "EG_transcript_2898", "EG_transcript_21535", "EG_transcript_20574", "EG_transcript_253", "EG_transcript_11544", "EG_transcript_6250", "EG_transcript_30017", "EG_transcript_3093", "EG_transcript_7320", "EG_transcript_50867", "EG_transcript_5289", "EG_transcript_9285", "EG_transcript_12637", "EG_transcript_37730", "EG_transcript_4177", "EG_transcript_17694", "EG_transcript_15128", "EG_transcript_29973", "EG_transcript_41717", "EG_transcript_17874", "EG_transcript_11626", "EG_transcript_260", "EG_transcript_13668", "EG_transcript_4738", "EG_transcript_15806", "EG_transcript_25802", "EG_transcript_8430", "EG_transcript_8430", "EG_transcript_8430", "EG_transcript_3260", "EG_transcript_54525", "EG_transcript_24888", "EG_transcript_38515", "EG_transcript_30675", "EG_transcript_37574", "EG_transcript_53416", "EG_transcript_1975", "EG_transcript_5491", "EG_transcript_18163", "EG_transcript_7825", "EG_transcript_1813", "EG_transcript_12219", "EG_transcript_25958", "EG_transcript_4296", "EG_transcript_10777", "EG_transcript_17880", "EG_transcript_12428", "EG_transcript_4312", "EG_transcript_17899", "EG_transcript_5801", "EG_transcript_22347", "EG_transcript_5443", "EG_transcript_6218", "EG_transcript_6315", "EG_transcript_22270", "EG_transcript_7318", "EG_transcript_12438", "EG_transcript_20497", "EG_transcript_15350", "EG_transcript_9396", "EG_transcript_18414", "EG_transcript_9730", "EG_transcript_22223", "EG_transcript_28843", "EG_transcript_9456", "EG_transcript_9184", "EG_transcript_13483", "EG_transcript_19314", "EG_transcript_10", "EG_transcript_21310", "EG_transcript_13992", "EG_transcript_7359", "EG_transcript_7896", "EG_transcript_298", "EG_transcript_4260", "EG_transcript_44081", "EG_transcript_34406", "EG_transcript_6269", "EG_transcript_18536", "EG_transcript_11871", "EG_transcript_4260", "EG_transcript_6510", "EG_transcript_17772", "EG_transcript_2406", "EG_transcript_10565", "EG_transcript_4382", "EG_transcript_35793", "EG_transcript_11", "EG_transcript_21", "EG_transcript_4390", "EG_transcript_9957", "EG_transcript_3246", "EG_transcript_22333", "EG_transcript_8430", "EG_transcript_13777", "EG_transcript_25414", "EG_transcript_21879", "EG_transcript_25617", "EG_transcript_6595", "EG_transcript_10242", "EG_transcript_21906", "EG_transcript_48180", "EG_transcript_36233", "EG_transcript_10372", "EG_transcript_8785", "EG_transcript_29409", "EG_transcript_19156", "EG_transcript_7630", "EG_transcript_18060", "EG_transcript_10541", "EG_transcript_17241", "EG_transcript_18103", "EG_transcript_8676", "EG_transcript_19552", "EG_transcript_30283", "EG_transcript_8117", "EG_transcript_14772", "EG_transcript_5340", "EG_transcript_41526", "EG_transcript_8180", "EG_transcript_8477", "EG_transcript_29203", "EG_transcript_1759", "EG_transcript_35745", "EG_transcript_7126", "EG_transcript_20283", "EG_transcript_23555", "EG_transcript_16062", "EG_transcript_3471", "EG_transcript_16257", "EG_transcript_13959", "EG_transcript_17531", "EG_transcript_3478", "EG_transcript_33117", "EG_transcript_19356", "EG_transcript_2238", "EG_transcript_42708", "EG_transcript_10201", "EG_transcript_15615", "EG_transcript_41901", "EG_transcript_7383", "EG_transcript_19909", "EG_transcript_4413", "EG_transcript_1264", "EG_transcript_4115", "EG_transcript_7672", "EG_transcript_9844", "EG_transcript_23810", "EG_transcript_12653", "EG_transcript_7474", "EG_transcript_13551", "EG_transcript_6158", "EG_transcript_5646", "EG_transcript_7383", "EG_transcript_5998", "EG_transcript_11375", "EG_transcript_19503", "EG_transcript_52704", "EG_transcript_19795", "EG_transcript_18", "EG_transcript_20144", "EG_transcript_16", "EG_transcript_24694", "EG_transcript_18914", "EG_transcript_14", "EG_transcript_9505", "EG_transcript_40935", "EG_transcript_20545", "EG_transcript_33577", "EG_transcript_27692", "EG_transcript_25501", "EG_transcript_13570", "EG_transcript_29149", "EG_transcript_13768", "EG_transcript_433", "EG_transcript_2978", "EG_transcript_18315", "EG_transcript_32395", "EG_transcript_12972", "EG_transcript_11180", "EG_transcript_59561", "EG_transcript_25604", "EG_transcript_18669", "EG_transcript_12684", "EG_transcript_39815", "EG_transcript_13256", "EG_transcript_22072", "EG_transcript_4666", "EG_transcript_23305", "EG_transcript_16257", "EG_transcript_29", "EG_transcript_18831", "EG_transcript_15548", "EG_transcript_18355", "EG_transcript_32956", "EG_transcript_25593", "EG_transcript_24567", "EG_transcript_17479", "EG_transcript_8095", "EG_transcript_6808", "EG_transcript_40714", "EG_transcript_11341", "EG_transcript_36862", "EG_transcript_31449", "EG_transcript_3873", "EG_transcript_16359", "EG_transcript_10562", "EG_transcript_17968", "EG_transcript_18464", "EG_transcript_7085", "EG_transcript_50540", "EG_transcript_6932", "EG_transcript_8912", "EG_transcript_13997", "EG_transcript_5315", "EG_transcript_8513", "EG_transcript_21", "EG_transcript_976", "EG_transcript_5584", "EG_transcript_19383", "EG_transcript_8552", "EG_transcript_27763", "EG_transcript_933", "EG_transcript_9931", "EG_transcript_30862", "EG_transcript_9915", "EG_transcript_18974", "EG_transcript_9483", "EG_transcript_42708", "EG_transcript_12885", "EG_transcript_18272", "EG_transcript_7094", "EG_transcript_15060", "EG_transcript_26048", "EG_transcript_14637", "EG_transcript_4309", "EG_transcript_9849", "EG_transcript_18070", "EG_transcript_3065", "EG_transcript_24767", "EG_transcript_32046", "EG_transcript_27781", "EG_transcript_8939", "EG_transcript_13", "EG_transcript_14515", "EG_transcript_39302", "EG_transcript_5897", "EG_transcript_11231", "EG_transcript_32996", "EG_transcript_4547", "EG_transcript_23844", "EG_transcript_31840", "EG_transcript_41871", "EG_transcript_4696", "EG_transcript_7801", "EG_transcript_4223", "EG_transcript_13577", "EG_transcript_8297", "EG_transcript_40617", "EG_transcript_13256", "EG_transcript_8672", "EG_transcript_15", "EG_transcript_19911", "EG_transcript_18", "EG_transcript_11996", "EG_transcript_2843", "EG_transcript_17851", "EG_transcript_62", "EG_transcript_10", "EG_transcript_45706", "EG_transcript_45706", "EG_transcript_30936", "EG_transcript_32001", "EG_transcript_22468", "EG_transcript_12", "EG_transcript_65", "EG_transcript_20434", "EG_transcript_20789", "EG_transcript_5237", "EG_transcript_2988", "EG_transcript_33225", "EG_transcript_2076", "EG_transcript_30796", "EG_transcript_15827", "EG_transcript_14", "EG_transcript_26376", "EG_transcript_724", "EG_transcript_15810", "EG_transcript_16", "EG_transcript_33134", "EG_transcript_5418", "EG_transcript_5418", "EG_transcript_11", "EG_transcript_33323", "EG_transcript_5645", "EG_transcript_15653", "EG_transcript_5475", "EG_transcript_19472", "EG_transcript_29012", "EG_transcript_17526", "EG_transcript_18691", "EG_transcript_11102", "EG_transcript_35867", "EG_transcript_4078", "EG_transcript_23444", "EG_transcript_14923", "EG_transcript_41155", "EG_transcript_22381", "EG_transcript_15742", "EG_transcript_30514", "EG_transcript_13613", "EG_transcript_18807", "EG_transcript_17654", "EG_transcript_1173", "EG_transcript_14855", "EG_transcript_26207", "EG_transcript_13", "EG_transcript_30891", "EG_transcript_28375", "EG_transcript_19208", "EG_transcript_25720", "EG_transcript_661", "EG_transcript_709", "EG_transcript_19311", "EG_transcript_9513", "EG_transcript_74", "EG_transcript_52320", "EG_transcript_14838", "EG_transcript_26395", "EG_transcript_21706", "EG_transcript_63", "EG_transcript_1772", "EG_transcript_30735", "EG_transcript_10179", "EG_transcript_24975", "EG_transcript_63", "EG_transcript_13881", "EG_transcript_19719", "EG_transcript_994", "EG_transcript_10609", "EG_transcript_23216", "EG_transcript_29710", "EG_transcript_11980", "EG_transcript_10370", "EG_transcript_9026", "EG_transcript_22409", "EG_transcript_51504", "EG_transcript_5071", "EG_transcript_2795", "EG_transcript_19750", "EG_transcript_17195", "EG_transcript_6492", "EG_transcript_20934", "EG_transcript_264", "EG_transcript_23080", "EG_transcript_14995", "EG_transcript_737", "EG_transcript_12763", "EG_transcript_8744", "EG_transcript_8850", "EG_transcript_3057", "EG_transcript_13068", "EG_transcript_3914", "EG_transcript_33", "EG_transcript_7147", "EG_transcript_28800", "EG_transcript_26412", "EG_transcript_74", "EG_transcript_2552", "EG_transcript_26284", "EG_transcript_11772", "EG_transcript_19677", "EG_transcript_39755", "EG_transcript_1851", "EG_transcript_17", "EG_transcript_14106", "EG_transcript_14106", "EG_transcript_12", "EG_transcript_39015", "EG_transcript_18209", "EG_transcript_7155", "EG_transcript_6511", "EG_transcript_1927", "EG_transcript_9", "EG_transcript_12195", "EG_transcript_524", "EG_transcript_23431", "EG_transcript_34", "EG_transcript_29", "EG_transcript_9069", "EG_transcript_21237", "EG_transcript_6513", "EG_transcript_20451", "EG_transcript_13807", "EG_transcript_47615", "EG_transcript_26203", "EG_transcript_26203", "EG_transcript_37499", "EG_transcript_106", "EG_transcript_8629", "EG_transcript_1615", "EG_transcript_29249", "EG_transcript_15055", "EG_transcript_52", "EG_transcript_4298", "EG_transcript_14625", "EG_transcript_33", "EG_transcript_3722", "EG_transcript_33121", "EG_transcript_26805", "EG_transcript_40131", "EG_transcript_52", "EG_transcript_2545", "EG_transcript_6401", "EG_transcript_26638", "EG_transcript_19809", "EG_transcript_9521", "EG_transcript_29902", "EG_transcript_2480", "EG_transcript_201", "EG_transcript_24541", "EG_transcript_22816", "EG_transcript_15072", "EG_transcript_18240", "EG_transcript_6189", "EG_transcript_34", "EG_transcript_7299", "EG_transcript_8750", "EG_transcript_150", "EG_transcript_7667", "EG_transcript_159", "EG_transcript_17416", "EG_transcript_1640", "EG_transcript_9635", "EG_transcript_19458", "EG_transcript_19184", "EG_transcript_6744", "EG_transcript_5202", "EG_transcript_5214", "EG_transcript_23831", "EG_transcript_7817", "EG_transcript_14371", "EG_transcript_3072", "EG_transcript_42748", "EG_transcript_10177", "EG_transcript_196", "EG_transcript_29941", "EG_transcript_10882", "EG_transcript_9967", "EG_transcript_26156", "EG_transcript_982", "EG_transcript_22169", "EG_transcript_22412", "EG_transcript_7177", "EG_transcript_1349", "EG_transcript_25719", "EG_transcript_26134", "EG_transcript_21174", "EG_transcript_11734", "EG_transcript_7177", "EG_transcript_3235", "EG_transcript_1067", "EG_transcript_21460", "EG_transcript_9431", "EG_transcript_18154", "EG_transcript_15833", "EG_transcript_17", "EG_transcript_364", "EG_transcript_37274", "EG_transcript_29972", "EG_transcript_31130", "EG_transcript_760", "EG_transcript_5230", "EG_transcript_10136", "EG_transcript_37022", "EG_transcript_13840", "EG_transcript_17771", "EG_transcript_364", "EG_transcript_15107", "EG_transcript_106", "EG_transcript_19005", "EG_transcript_372", "EG_transcript_3497", "EG_transcript_31367", "EG_transcript_3096", "EG_transcript_29479", "EG_transcript_12231", "EG_transcript_372", "EG_transcript_526", "EG_transcript_14405", "EG_transcript_10463", "EG_transcript_35367", "EG_transcript_1740", "EG_transcript_12216", "EG_transcript_10194", "EG_transcript_39280", "EG_transcript_162", "EG_transcript_21390", "EG_transcript_16113", "EG_transcript_13588", "EG_transcript_13915", "EG_transcript_56241", "EG_transcript_13555", "EG_transcript_6098", "EG_transcript_2814", "EG_transcript_20669", "EG_transcript_17", "EG_transcript_2593", "EG_transcript_9247", "EG_transcript_10550", "EG_transcript_27003", "EG_transcript_29299", "EG_transcript_2731", "EG_transcript_1661", "EG_transcript_23125", "EG_transcript_21060"]
small_q = ["Tb927.10.12840", "Tb927.10.12840", "Tb927.10.14820", "Tb927.10.14820", "Tb927.10.14820", "Tb927.10.14830", "Tb927.10.14830", "Tb927.10.14830", "Tb927.10.14840", "Tb927.10.14840", "Tb927.10.14840", "Tb927.10.15790", "Tb927.10.2230", "Tb927.10.2350", "Tb927.10.2350", "Tb927.10.2350", "Tb927.10.2560", "Tb927.10.2960", "Tb927.10.4130", "Tb927.10.4280", "Tb927.10.4910", "Tb927.10.6400", "Tb927.10.6400", "Tb927.10.6510", "Tb927.10.6510", "Tb927.10.7570", "Tb927.10.7570", "Tb927.11.1310", "Tb927.11.13140", "Tb927.11.13580", "Tb927.11.13750", "Tb927.11.14440", "Tb927.11.15040", "Tb927.11.15040", "Tb927.11.15240", "Tb927.11.15240", "Tb927.11.15550", "Tb927.11.180", "Tb927.11.270", "Tb927.11.3750", "Tb927.11.3980", "Tb927.11.530", "Tb927.11.6140", "Tb927.11.8380", "Tb927.11.9900", "Tb927.2.1560", "Tb927.2.1560", "Tb927.2.1680", "Tb927.2.1680", "Tb927.2.2510", "Tb927.2.2520", "Tb927.3.1000", "Tb927.3.2300", "Tb927.4.1660", "Tb927.5.1060", "Tb927.6.2420", "Tb927.6.4980", "Tb927.7.2700", "Tb927.7.3520", "Tb927.7.4140", "Tb927.7.7230", "Tb927.7.7420", "Tb927.7.7430", "Tb927.8.1890", "Tb927.8.4330", "Tb927.8.4330", "Tb927.8.4330", "Tb927.8.4610", "Tb927.8.4810", "Tb927.8.5120", "Tb927.8.5120", "Tb927.8.5810", "Tb927.8.5810", "Tb927.8.6580", "Tb927.8.890", "Tb927.8.890", "Tb927.9.11040", "Tb927.9.11040", "Tb927.9.14160", "Tb927.9.15240", "Tb927.9.4520"]
big_q = ["Tb927.8.2050", "Tb927.3.4380", "Tb927.9.10770", "Tb927.11.1270", "Tb927.10.8900", "Tb927.7.4810", "Tb927.10.5840", "Tb927.11.3750", "Tb927.9.15380", "Tb927.1.2430", "Tb927.10.12500", "Tb927.11.9560", "Tb927.4.2310", "Tb927.11.2650", "Tb927.2.5210", "Tb927.11.3730", "Tb927.8.2030", "Tb927.8.6750", "Tb927.1.710", "Tb927.2.3080", "Tb927.11.900", "Tb927.10.14780", "Tb927.3.2310", "Tb927.9.4310", "Tb927.10.12510", "Tb927.10.4560", "Tb927.8.4330", "Tb927.10.12500", "Tb927.10.12500", "Tb927.3.930", "Tb927.5.1520", "Tb927.11.13440", "Tb927.10.2560", "Tb927.11.7460", "Tb927.8.7100", "Tb927.11.6140", "Tb927.11.9720", "Tb927.7.3410", "Tb927.11.16730", "Tb927.10.12630", "Tb927.9.11040", "Tb927.11.13750", "Tb927.11.1430", "Tb927.11.11780", "Tb927.10.12500", "Tb927.10.16120", "Tb927.10.12510", "Tb927.10.180", "Tb927.7.1900", "Tb927.8.1240", "Tb927.1.3110", "Tb927.11.5520", "Tb927.7.1910", "Tb927.2.4240", "Tb927.11.16480", "Tb927.7.210", "Tb927.9.12550", "Tb927.9.12570", "Tb927.9.12630", "Tb927.9.12590", "Tb927.9.12610", "Tb927.10.10390", "Tb927.11.13650", "Tb927.10.4910", "Tb927.10.14000", "Tb927.10.14000", "Tb927.10.6400", "Tb927.10.6510", "Tb927.10.12510", "Tb927.9.8950", "Tb927.3.4290", "Tb927.6.2490", "Tb927.9.12650", "Tb927.5.900", "Tb927.3.4190", "Tb927.3.930", "Tb927.11.5050", "Tb927.3.3900", "Tb927.10.11930", "Tb927.8.1860", "Tb927.7.3500", "Tb927.11.6230", "Tb927.10.2100", "Tb927.10.2530", "Tb927.3.4650", "Tb927.8.1160", "Tb927.10.4560", "Tb927.8.6060", "Tb927.10.1070", "Tb927.8.6110", "Tb927.10.12000", "Tb927.10.14280", "Tb927.7.4480", "Tb927.6.4280", "Tb927.5.1130", "Tb927.11.9670", "Tb927.11.2360", "Tb927.7.6930", "Tb927.9.4680", "Tb927.11.5290", "Tb927.6.4280", "Tb927.10.4990", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.2.4130", "Tb927.9.5690", "Tb927.11.14440", "Tb927.11.7170", "Tb927.11.6280", "Tb927.4.4680", "Tb927.10.6050", "Tb927.7.4770", "Tb927.6.3630", "Tb927.7.4160", "Tb927.10.11390", "Tb927.10.8630", "Tb927.10.14780", "Tb927.1.3950", "Tb927.11.14120", "Tb927.9.7770", "Tb927.3.1790", "Tb927.11.15860", "Tb927.11.15840", "Tb927.9.5900", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.10.4000", "Tb927.10.4040", "Tb927.3.4850", "Tb927.10.4760", "Tb927.10.14150", "Tb927.11.15750", "Tb927.3.4290", "Tb927.11.5440", "Tb927.8.5810", "Tb927.10.12970", "Tb927.11.8090", "Tb927.10.7410", "Tb927.6.4540", "Tb927.9.7110", "Tb927.10.6190", "Tb927.10.6190", "Tb927.10.6400", "Tb927.10.6510", "Tb11.02.5380", "Tb927.10.2890", "Tb927.9.12730", "Tb927.6.2170", "Tb927.3.3330", "Tb927.1.2330", "Tb927.8.2030", "Tb927.9.15090", "Tb927.10.660", "Tb927.2.3030", "Tb927.11.4700", "Tb927.7.5720", "Tb927.9.6310", "Tb927.7.6260", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.2.3030", "Tb927.10.12700", "Tb927.4.480", "Tb927.8.4970", "Tb927.8.8020", "Tb927.9.1780", "Tb927.11.11330", "Tb927.11.1540", "Tb927.11.7380", "Tb927.2.4110", "Tb927.10.6910", "Tb927.8.2110", "Tb927.8.4970", "Tb927.10.13480", "Tb927.8.6580", "Tb927.11.7460", "Tb927.9.4950", "Tb927.11.11820", "Tb927.8.1160", "Tb927.7.1300", "Tb927.11.6210", "Tb927.10.3120", "Tb927.4.1660", "Tb927.7.1130", "Tb927.1.1770", "Tb927.10.9190", "Tb927.4.2260", "Tb927.8.3690", "Tb927.6.1570", "Tb927.11.2000", "Tb927.11.7460", "Tb927.11.11330", "Tb927.11.190", "Tb927.4.2700", "Tb927.11.7460", "Tb927.11.270", "Tb927.9.9740", "Tb927.11.11330", "Tb927.4.5040", "Tb927.7.7420", "Tb927.7.7430", "Tb927.1.1000", "Tb927.7.5480", "Tb927.10.13360", "Tb927.10.8730", "Tb927.4.3660", "Tb927.5.3400", "Tb927.8.6240", "Tb927.10.5620", "Tb927.9.4230", "Tb927.6.4980", "Tb927.7.2700", "Tb927.10.10310", "Tb927.8.5120", "Tb927.6.5080", "Tb927.2.4590", "Tb927.2.4610", "Tb927.11.15530", "Tb927.3.2230", "Tb927.11.8380", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.9.3350", "Tb927.10.14840", "Tb927.10.14830", "Tb927.10.14820", "Tb927.8.4580", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.3.3580", "Tb927.11.1670", "Tb927.5.800", "Tb927.11.9860", "Tb927.7.5540", "Tb927.11.11290", "Tb927.11.11950", "Tb927.10.5560", "Tb927.5.2380", "Tb927.6.2600", "Tb927.9.6230", "Tb927.6.3050", "Tb927.8.1890", "Tb927.3.2980", "Tb927.3.1380", "Tb927.11.10430", "Tb927.11.11730", "Tb927.11.5680", "Tb927.2.2970", "Tb927.10.15490", "Tb927.9.5890", "Tb927.11.15750", "Tb927.10.3080", "Tb927.9.9740", "Tb927.3.680", "Tb927.9.9740", "Tb927.6.2170", "Tb927.11.15750", "Tb927.10.5330", "Tb927.8.2670", "Tb927.9.9740", "Tb927.9.2560", "Tb927.11.15040", "Tb927.6.2540", "Tb927.10.8230", "Tb927.2.1680", "Tb927.2.1560", "Tb927.7.1290", "Tb927.9.6310", "Tb927.10.4310", "Tb10.v4.0045", "Tb927.3.3270", "Tb927.8.7980", "Tb927.4.4380", "Tb927.4.1660", "Tb927.11.15040", "Tb927.7.5770", "Tb927.11.14360", "Tb927.1.1770", "Tb927.11.13580", "Tb927.10.12890", "Tb927.8.3530", "Tb927.9.6170", "Tb927.2.3030", "Tb927.8.5460", "Tb927.10.6730", "Tb927.10.6030", "Tb927.9.3990", "Tb927.11.13220", "Tb927.5.1210", "Tb927.2.6070", "Tb927.10.4330", "Tb927.8.1160", "Tb927.8.1160", "Tb927.7.3750", "Tb927.11.14350", "Tb927.11.3250", "Tb927.8.1720", "Tb927.10.15010", "Tb927.10.13040", "Tb927.5.1470", "Tb927.11.16130", "Tb927.11.16130", "Tb927.11.5450", "Tb927.2.3030", "Tb927.11.1370", "Tb927.9.6040", "Tb927.10.6060", "Tb927.7.3940", "Tb927.11.3250", "Tb927.9.3170", "Tb927.10.6610", "Tb927.11.4780", "Tb927.6.700", "Tb927.11.7460", "Tb927.7.4950", "Tb927.8.1160", "Tb927.5.450", "Tb927.11.15750", "Tb927.11.15550", "Tb927.6.2060", "Tb927.11.16740", "Tb927.11.12590", "Tb927.8.6060", "Tb927.3.3630", "Tb927.7.1640", "Tb927.11.15650", "Tb927.11.15020", "Tb927.10.8050", "Tb927.9.1960", "Tb927.10.3940", "Tb927.11.2650", "Tb927.11.2650", "Tb927.9.4210", "Tb927.11.10510", "Tb927.11.11330", "Tb927.10.7700", "Tb927.2.3030", "Tb927.10.13430", "Tb927.8.2500", "Tb927.11.15240", "Tb927.11.11540", "Tb927.9.10310", "Tb927.8.7120", "Tb927.9.11000", "Tb927.4.1280", "Tb927.8.2240", "Tb927.10.3370", "Tb927.11.14480", "Tb927.9.3370", "Tb927.11.3980", "Tb927.8.920", "Tb927.4.3300", "Tb927.9.9710", "Tb927.6.1550", "Tb927.7.1130", "Tb927.11.13020", "Tb927.9.10080", "Tb927.7.7410", "Tb927.8.2160", "Tb927.8.4640", "Tb927.7.4940", "Tb927.3.2840", "Tb927.11.10240", "Tb927.6.3740", "Tb927.6.3750", "Tb927.6.3800", "Tb927.11.16930", "Tb927.10.15790", "Tb927.8.4810", "Tb927.9.3370", "Tb927.9.9840", "Tb927.5.4170", "Tb927.5.4170", "Tb927.10.7700", "Tb927.6.4210", "Tb927.5.3220", "Tb927.10.2440", "Tb927.5.3010", "Tb927.8.2400", "Tb927.11.3590", "Tb927.10.3260", "Tb927.7.220", "Tb927.11.8990", "Tb927.10.4940", "Tb927.10.5010", "Tb927.10.8110", "Tb927.10.2960", "Tb927.3.1120", "Tb927.3.1590", "Tb927.9.4190", "Tb927.6.3510", "Tb927.7.4710", "Tb927.8.8300", "Tb927.8.3580", "Tb927.5.2590", "Tb927.10.12840", "Tb927.10.9820", "Tb927.9.12120", "Tb927.4.3680", "Tb927.2.5180", "Tb927.6.4080", "Tb927.5.960", "Tb927.11.13180", "Tb927.8.3330", "Tb927.2.2430", "Tb927.3.930", "Tb927.9.11720", "Tb927.8.2400", "Tb927.8.2400", "Tb927.10.14090", "Tb927.11.17000", "Tb927.11.7460", "Tb927.7.2820", "Tb927.7.2820", "Tb927.8.2850", "Tb927.10.13850", "Tb927.7.7360", "Tb927.11.11330", "Tb927.11.6630", "Tb927.9.2320", "Tb927.11.8870", "Tb927.7.5120", "Tb927.7.3280", "Tb927.10.1100", "Tb927.3.930", "Tb927.11.3250", "Tb927.8.4000", "Tb927.3.4260", "Tb927.6.700", "Tb927.7.4170", "Tb927.11.11290", "Tb927.11.15950", "Tb927.11.15130", "Tb927.3.1730", "Tb927.11.7540", "Tb927.10.9900", "Tb927.10.8490", "Tb927.9.10580", "Tb927.10.10460", "Tb927.4.2740", "Tb927.8.2400", "Tb927.11.3220", "Tb927.7.890", "Tb927.11.10140", "Tb927.8.2920", "Tb927.2.6070", "Tb927.8.1020", "Tb927.11.10020", "Tb927.10.12240", "Tb927.5.3420", "Tb927.11.9390", "Tb927.10.10110", "Tb927.1.1000", "Tb927.8.8360", "Tb927.8.560", "Tb927.2.4700", "Tb927.7.2710", "Tb927.8.7730", "Tb927.7.5790", "Tb927.10.7700", "Tb927.5.3350", "Tb927.7.2410", "Tb927.7.3430", "Tb927.10.12960", "Tb927.10.15420", "Tb927.6.1510", "Tb927.8.7530", "Tb927.11.5090", "Tb927.4.5010", "Tb927.7.190", "Tb927.11.15230", "Tb927.5.1310", "Tb927.9.9860", "Tb927.9.5750", "Tb927.10.13120", "Tb927.11.6000", "Tb927.1.1770", "Tb927.11.7460", "Tb927.8.2400", "Tb927.11.6460", "Tb927.11.5970", "Tb927.11.540", "Tb927.10.1260", "Tb927.5.2780", "Tb927.7.4420", "Tb927.9.8410", "Tb927.10.8450", "Tb927.6.2790", "Tb927.11.3900", "Tb927.10.3210", "Tb927.11.11330", "Tb927.11.2920", "Tb927.10.2000", "Tb927.11.530", "Tb927.7.1040", "Tb927.7.2190", "Tb927.3.930", "Tb927.7.5680", "Tb927.3.930", "Tb927.6.2490", "Tb927.8.6390", "Tb927.3.930", "Tb927.11.3010", "Tb927.10.11300", "Tb927.8.7170", "Tb927.11.180", "Tb927.7.1720", "Tb927.6.2420", "Tb927.9.4620", "Tb927.10.11300", "Tb927.4.2480", "Tb927.10.16000", "Tb927.11.6660", "Tb927.10.7090", "Tb927.11.1320", "Tb927.11.14980", "Tb927.5.4330", "Tb927.10.1580", "Tb927.10.11900", "Tb927.10.10420", "Tb927.4.1540", "Tb927.6.710", "Tb927.7.610", "Tb927.7.7230", "Tb927.10.10610", "Tb927.10.710", "Tb927.4.4910", "Tb927.3.930", "Tb927.11.15820", "Tb927.6.1140", "Tb927.6.4180", "Tb927.6.4990", "Tb927.10.9270", "Tb927.9.15240", "Tb927.1.1770", "Tb927.10.7620", "Tb927.11.4490", "Tb927.8.890", "Tb927.9.11540", "Tb927.5.1030", "Tb927.4.1070", "Tb927.10.1500", "Tb927.10.12880", "Tb927.1.2340", "Tb927.7.3520", "Tb927.11.8800", "Tb927.7.1030", "Tb927.11.6300", "Tb927.8.3060", "Tb927.3.1380", "Tb927.10.4280", "Tb927.2.3030", "Tb927.10.1050", "Tb927.3.930", "Tb927.8.650", "Tb927.10.7570", "Tb927.11.480", "Tb927.11.6150", "Tb927.2.5210", "Tb927.5.320", "Tb927.8.7600", "Tb927.11.7040", "Tb927.6.2010", "Tb927.3.5610", "Tb927.10.10160", "Tb927.8.1990", "Tb927.10.12540", "Tb927.9.2670", "Tb927.11.12350", "Tb927.6.1010", "Tb927.8.6390", "Tb927.9.14070", "Tb927.4.4210", "Tb927.10.13860", "Tb927.6.1250", "Tb927.11.1680", "Tb927.3.4990", "Tb927.3.800", "Tb927.10.12930", "Tb927.4.5390", "Tb927.3.930", "Tb927.6.4560", "Tb927.11.13620", "Tb927.7.4730", "Tb927.11.3700", "Tb927.2.4090", "Tb927.11.530", "Tb927.9.14160", "Tb927.10.14600", "Tb927.5.1710", "Tb927.7.6460", "Tb927.6.4320", "Tb927.4.2950", "Tb927.10.9860", "Tb927.11.4480", "Tb927.9.3590", "Tb927.7.600", "Tb927.9.7190", "Tb927.11.3250", "Tb927.3.4820", "Tb927.11.3250", "Tb927.10.11940", "Tb927.6.4760", "Tb927.9.12100", "Tb927.4.420", "Tb927.11.3250", "Tb927.11.5280", "Tb927.10.1570", "Tb927.10.6120", "Tb927.3.1000", "Tb927.3.2150", "Tb927.11.3250", "Tb927.8.7100", "Tb927.9.8680", "Tb927.2.3460", "Tb927.10.10440", "Tb927.7.4180", "Tb927.11.8910", "Tb927.5.320", "Tb927.8.4610", "Tb927.11.15150", "Tb927.11.3250", "Tb927.7.3960", "Tb927.8.6970", "Tb927.11.9900", "Tb927.11.3250", "Tb927.9.15360", "Tb927.6.2350", "Tb927.6.2290", "Tb927.11.3250", "Tb927.9.5960", "Tb927.3.3560", "Tb927.8.5450", "Tb927.7.740", "Tb927.7.6350", "Tb927.9.7110", "Tb927.9.7110", "Tb927.1.740", "Tb927.8.5140", "Tb927.4.4980", "Tb927.9.15460", "Tb927.7.4440", "Tb927.11.11440", "Tb927.3.860", "Tb927.11.7212", "Tb927.5.520", "Tb927.9.7110", "Tb927.8.6410", "Tb927.7.6800", "Tb927.11.9360", "Tb927.7.4460", "Tb927.9.2450", "Tb927.6.2510", "Tb927.11.3250", "Tb927.11.6640", "Tb927.5.1940", "Tb927.8.5540", "Tb927.8.3380", "Tb927.9.10640", "Tb927.8.7290", "Tb927.3.5360", "Tb927.9.10010", "Tb927.3.930", "Tb927.7.4140", "Tb927.8.5140", "Tb927.7.4910", "Tb927.5.1060", "Tb927.3.930", "Tb927.7.1440", "Tb927.11.13230", "Tb927.8.3130", "Tb927.10.13300", "Tb927.11.3250", "Tb927.11.5240", "Tb927.10.13620", "Tb927.8.2630", "Tb927.4.700", "Tb927.10.2190", "Tb927.11.14080", "Tb927.4.2450", "Tb927.1.1330", "Tb927.8.8300", "Tb927.10.8030", "Tb927.10.5610", "Tb927.11.9310", "Tb927.3.5130", "Tb927.2.3460", "Tb927.8.2970", "Tb927.5.4070", "Tb927.10.14740", "Tb927.10.14010", "Tb927.8.2400", "Tb927.10.3690", "Tb927.6.1680", "Tb927.7.910", "Tb927.8.3830", "Tb927.10.13510", "Tb927.11.13760", "Tb927.11.3900", "Tb927.10.8020", "Tb927.3.930", "Tb927.8.6640", "Tb927.11.7900", "Tb927.5.1510", "Tb927.11.3250", "Tb927.8.700", "Tb927.6.1840", "Tb927.8.7730", "Tb927.11.9820", "Tb927.11.16430", "Tb927.1.4050", "Tb927.3.930", "Tb927.9.6090", "Tb927.9.6100", "Tb927.3.930", "Tb927.11.1310", "Tb927.8.6820", "Tb927.7.5820", "Tb927.7.900", "Tb927.6.3350", "Tb927.9.11900", "Tb927.11.15150", "Tb927.9.3100", "Tb927.10.6850", "Tb927.11.3250", "Tb927.11.3250", "Tb927.10.12980", "Tb927.6.4930", "Tb927.9.5590", "Tb927.9.9150", "Tb927.10.6300", "Tb927.10.4880", "Tb927.2.4740", "Tb927.2.4890", "Tb927.4.1130", "Tb927.3.930", "Tb927.9.9150", "Tb927.1.4050", "Tb927.11.5180", "Tb927.9.13580", "Tb927.3.930", "Tb927.10.9080", "Tb927.3.1890", "Tb927.11.3250", "Tb927.8.7290", "Tb927.11.530", "Tb927.9.14200", "Tb927.3.2300", "Tb927.11.3250", "Tb927.11.5780", "Tb927.9.2620", "Tb927.1.4100", "Tb927.9.2450", "Tb927.7.1550", "Tb927.11.5820", "Tb927.11.14440", "Tb927.8.4540", "Tb927.10.4130", "Tb927.9.2050", "Tb927.8.4540", "Tb927.11.4200", "Tb927.10.5590", "Tb927.3.930", "Tb927.4.1470", "Tb927.7.4310", "Tb927.11.3250", "Tb927.7.940", "Tb927.1.880", "Tb927.11.2910", "Tb927.10.9970", "Tb927.11.5740", "Tb927.11.530", "Tb927.3.1840", "Tb927.11.3940", "Tb927.1.2990", "Tb927.4.3950", "Tb927.8.2540", "Tb927.1.2120", "Tb927.10.16150", "Tb927.11.12220", "Tb927.7.4550", "Tb927.11.7960", "Tb927.3.4630", "Tb927.5.1550", "Tb927.11.6980", "Tb927.5.3360", "Tb927.8.2540", "Tb927.9.4190", "Tb927.7.2620", "Tb927.9.11280", "Tb927.10.3650", "Tb927.4.590", "Tb927.9.4520", "Tb927.8.5860", "Tb927.3.2130", "Tb927.11.11460", "Tb927.5.930", "Tb927.2.3180", "Tb927.4.1500", "Tb927.11.10960", "Tb927.11.5740", "Tb927.10.2230", "Tb927.8.6420", "Tb927.11.3250", "Tb927.3.930", "Tb927.1.1200", "Tb927.8.630", "Tb927.11.13290", "Tb927.8.4950", "Tb927.8.7290", "Tb927.11.5740", "Tb927.9.12480", "Tb927.10.7520", "Tb927.10.13830", "Tb927.11.3250", "Tb927.10.6090", "Tb927.11.3250", "Tb927.6.2480", "Tb927.3.930", "Tb927.9.9630", "Tb927.10.680", "Tb927.11.5500", "Tb927.9.6060", "Tb927.11.4320", "Tb927.11.3250", "Tb927.3.930", "Tb927.8.2650", "Tb927.9.8160", "Tb927.10.2350", "Tb927.10.13740", "Tb927.10.12810", "Tb927.11.6040", "Tb927.5.1530", "Tb927.8.7290", "Tb927.11.3690", "Tb927.3.950", "Tb927.11.16940", "Tb927.9.10400", "Tb927.10.2350", "Tb927.11.15470", "Tb927.8.7290", "Tb927.11.10150", "Tb927.11.14730", "Tb927.3.930", "Tb927.10.7310", "Tb927.8.1060", "Tb927.5.2790", "Tb927.4.2450", "Tb927.8.4010", "Tb927.9.6710", "Tb927.7.3950", "Tb927.4.2450", "Tb927.11.11490"]

# print len(set(big_q))
results = []
for word in small_q:
	if word not in big_q:
		results.append(word)

for i in set(results):
	print i