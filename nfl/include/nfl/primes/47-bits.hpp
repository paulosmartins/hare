const unsigned int kMaxNbModuli = 37;

const value_type P[37] = { 140737471578113ULL, 140737454800897ULL, 140737414955009ULL, 140737398177793ULL, 140737383497729ULL, 140737352040449ULL, 140737314291713ULL, 140737282834433ULL, 140737251377153ULL, 140737175879681ULL, 140737171685377ULL, 140737121353729ULL, 140737089896449ULL, 140737083604993ULL, 140737031176193ULL, 140737014398977ULL, 140736982941697ULL, 140736976650241ULL, 140736943095809ULL, 140736901152769ULL, 140736867598337ULL, 140736848723969ULL, 140736842432513ULL, 140736825655297ULL, 140736762740737ULL, 140736724992001ULL, 140736716603393ULL, 140736706117633ULL, 140736653688833ULL, 140736628523009ULL, 140736599162881ULL, 140736555122689ULL, 140736548831233ULL, 140736471236609ULL, 140736458653697ULL, 140736439779329ULL, 140736416710657ULL};

// The associated lower word of their Newton quotients
const value_type Pn[37] = { 288230393331580927ULL, 576460872562532352ULL, 1261008536151061976ULL, 1549239247310511768ULL, 1801441175946583616ULL, 2341874057334555658ULL, 2990393833902384733ULL, 3530827246794317898ULL, 4071260901279219867ULL, 5368302657743266722ULL, 5440360573904820853ULL, 6305055902853290727ULL, 6845490797517982810ULL, 6953577805442191532ULL, 7854303247290268236ULL, 8142535530416918902ULL, 8682971246501744153ULL, 8791058418710045605ULL, 9367523500364889119ULL, 10088105238985016027ULL, 10664570939122615186ULL, 10988833016247555726ULL, 11096920394616817199ULL, 11385153450846808739ULL, 12466028023750963179ULL, 13114553231357015595ULL, 13258669991403759752ULL, 13438815965621778793ULL, 14339546239371892077ULL, 14771897009146826342ULL, 15276306435975515914ULL, 16032920970826191208ULL, 16141008800174609087ULL, 17474092823390361769ULL, 17690268749454496677ULL, 18014532711029866263ULL, 18410855448847339957ULL};

static constexpr unsigned int kModulusBitsize = 47;
const unsigned int kModulusRepresentationBitsize = 64;

// A primitive 2*kMaxPolyDegree root of unity for each one of the moduli
const value_type primitive_roots[37] = { 13082280195693ULL, 132249819143984ULL, 7657919137775ULL, 99266397934387ULL, 49450664300316ULL, 78078923262139ULL, 68585996457375ULL, 50416215073441ULL, 5402425042970ULL, 14680651199950ULL, 97250969347701ULL, 49040658083414ULL, 51917152515975ULL, 6193354483468ULL, 129626674339570ULL, 120208465818261ULL, 23910876808596ULL, 86757462093809ULL, 101418216873778ULL, 17623779479361ULL, 95047519205795ULL, 92152946364860ULL, 114914218281652ULL, 103402610267206ULL, 43553807091687ULL, 110300032686878ULL, 28016114902463ULL, 84138042815162ULL, 98416262504601ULL, 98964555607215ULL, 3411452020160ULL, 138349894563023ULL, 24414563367985ULL, 22694163034826ULL, 10160452556316ULL, 21147180895594ULL, 103757640955800ULL };

// Inverses of kMaxPolyDegree (for the other degrees it can be derived easily)
// for the different moduli
const value_type invkMaxPolyDegree[37] = { 140737337360401ULL, 140737320583201ULL, 140737280737351ULL, 140737263960151ULL, 140737249280101ULL, 140737217822851ULL, 140737180074151ULL, 140737148616901ULL, 140737117159651ULL, 140737041662251ULL, 140737037467951ULL, 140736987136351ULL, 140736955679101ULL, 140736949387651ULL, 140736896958901ULL, 140736880181701ULL, 140736848724451ULL, 140736842433001ULL, 140736808878601ULL, 140736766935601ULL, 140736733381201ULL, 140736714506851ULL, 140736708215401ULL, 140736691438201ULL, 140736628523701ULL, 140736590775001ULL, 140736582386401ULL, 140736571900651ULL, 140736519471901ULL, 140736494306101ULL, 140736464946001ULL, 140736420905851ULL, 140736414614401ULL, 140736337019851ULL, 140736324436951ULL, 140736305562601ULL, 140736282493951ULL };

// Polynomial related data
const unsigned int kMaxPolyDegree = 1048576;