# from TCutility.results.result import Result
import periodictable as pt

# class Radii:
#     def __init__(self, element, symbol, number, **properties):
#         self.element = element
#         self.symbol = symbol
#         self.number = number
#         for key, value in properties.items():
#             setattr(self, key, value)


element_order = ["hydrogen","helium","lithium","beryllium","boron","carbon","nitrogen","oxygen","fluorine","neon","sodium","magnesium","aluminium","silicon","phosphorus","sulfur","chlorine","argon","potassium","calcium","scandium","titanium","vanadium","chromium","manganese","iron","cobalt","nickel","copper","zinc","gallium","germanium","arsenic","selenium","bromine","krypton","rubidium","strontium","yttrium","zirconium","niobium","molybdenum","technetium","ruthenium","rhodium","palladium","silver","cadmium","indium","tin","antimony","tellurium","iodine","xenon","caesium","barium","lanthanum","cerium","praseodymium","neodymium","promethium","samarium","europium","gadolinium","terbium","dysprosium","holmium","erbium","thulium","ytterbium","lutetium","hafnium","tantalum","tungsten","rhenium","osmium","iridium","platinum","gold","mercury","thallium","lead","bismuth","polonium","astatine","radon","francium","radium","actinium","thorium","protactinium","uranium","neptunium","plutonium","americium","curium","berkelium","californium","einsteinium","fermium","mendelevium","nobelium","lawrencium","rutherfordium","dubnium","seaborgium","bohrium","hassium","meitnerium","darmstadtium","roentgenium","copernicium","nihonium","flerovium","moscovium","livermorium","tennessine","oganesson",]
symbol_order = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",]
colors = (0, 0, 0, 0), (1, 255, 255, 255), (2, 217, 255, 255), (3, 204, 128, 255), (4, 194, 255, 0), (5, 255, 181, 181), (6, 144, 144, 144), (7, 48, 80, 248), (8, 255, 13, 13), (9, 144, 224, 80), (10, 179, 227, 245), (11, 171, 92, 242), (12, 138, 255, 0), (13, 191, 166, 166), (14, 240, 200, 160), (15, 255, 128, 0), (16, 255, 255, 48), (17, 31, 240, 31), (18, 128, 209, 227), (19, 143, 64, 212), (20, 61, 255, 0), (21, 230, 230, 230), (22, 191, 194, 199), (23, 166, 166, 171), (24, 138, 153, 199), (25, 156, 122, 199), (26, 224, 102, 51), (27, 240, 144, 160), (28, 80, 208, 80), (29, 200, 128, 51), (30, 125, 128, 176), (31, 194, 143, 143), (32, 102, 143, 143), (33, 189, 128, 227), (34, 255, 161, 0), (35, 166, 41, 41), (36, 92, 184, 209), (37, 112, 46, 176), (38, 0, 255, 0), (39, 148, 255, 255), (40, 148, 224, 224), (41, 115, 194, 201), (42, 84, 181, 181), (43, 59, 158, 158), (44, 36, 143, 143), (45, 10, 125, 140), (46, 0, 105, 133), (47, 192, 192, 192), (48, 255, 217, 143), (49, 166, 117, 115), (50, 102, 128, 128), (51, 158, 99, 181), (52, 212, 122, 0), (53, 148, 0, 148), (54, 66, 158, 176), (55, 87, 23, 143), (56, 0, 201, 0), (57, 112, 212, 255), (58, 255, 255, 199), (59, 217, 255, 199), (60, 199, 255, 199), (61, 163, 255, 199), (62, 143, 255, 199), (63, 97, 255, 199), (64, 69, 255, 199), (65, 48, 255, 199), (66, 31, 255, 199), (67, 0, 255, 156), (68, 0, 230, 117), (69, 0, 212, 82), (70, 0, 191, 56), (71, 0, 171, 36), (72, 77, 194, 255), (73, 77, 166, 255), (74, 33, 148, 214), (75, 38, 125, 171), (76, 38, 102, 150), (77, 23, 84, 135), (78, 208, 208, 224), (79, 255, 209, 35), (80, 184, 184, 208), (81, 166, 84, 77), (82, 87, 89, 97), (83, 158, 79, 181), (84, 171, 92, 0), (85, 117, 79, 69), (86, 66, 130, 150), (87, 66, 0, 102), (88, 0, 125, 0), (89, 112, 171, 250), (90, 0, 186, 255), (91, 0, 161, 255), (92, 0, 143, 255), (93, 0, 128, 255), (94, 0, 107, 255), (95, 84, 92, 242), (96, 120, 92, 227), (97, 138, 79, 227), (98, 161, 54, 212), (99, 179, 31, 212), (100, 179, 31, 186), (101, 179, 13, 166), (102, 189, 13, 135), (103, 199, 0, 102), (104, 204, 0, 89), (105, 209, 0, 79), (106, 217, 0, 69), (107, 224, 0, 56), (108, 230, 0, 46), (109, 235, 0, 38), 


def parse_element(val):
    if isinstance(val, int):
        return val
    if val.lower() in element_order:
        return element_order.index(val.lower()) + 1
    if val in symbol_order:
        return symbol_order.index(val) + 1


def radius(element):
    num = parse_element(element)
    return pt.elements[num].covalent_radius


def color(element):
    num = parse_element(element)
    return colors[num][1:]


if __name__ == '__main__':
    print(radius(1))
