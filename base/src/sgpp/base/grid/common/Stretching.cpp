// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/Stretching.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

namespace sgpp {
namespace base {

static int leftIdx[2047] = {
  -2, -2, 0, -2, 1, 0, 2, -2, 3, 1, 4, 0, 5, 2, 6, -2, 7, 3, 8, 1, 9,
  4, 10, 0, 11, 5, 12, 2, 13, 6, 14, -2,
  15, 7, 16, 3, 17, 8, 18, 1, 19, 9, 20, 4, 21, 10, 22, 0, 23, 11, 24,
  5, 25, 12, 26, 2, 27, 13, 28, 6, 29, 14, 30, -2, 31, 15, 32,
  7, 33, 16, 34, 3, 35, 17, 36, 8, 37, 18, 38, 1, 39, 19, 40, 9, 41, 20,
  42, 4, 43, 21, 44, 10, 45, 22, 46, 0, 47, 23, 48, 11, 49,
  24, 50, 5, 51, 25, 52, 12, 53, 26, 54, 2, 55, 27, 56, 13, 57, 28, 58,
  6, 59, 29, 60, 14, 61, 30, 62, -2, 63, 31, 64, 15, 65, 32,
  66, 7, 67, 33, 68, 16, 69, 34, 70, 3, 71, 35, 72, 17, 73, 36, 74, 8,
  75, 37, 76, 18, 77, 38, 78, 1, 79, 39, 80, 19, 81, 40, 82, 9,
  83, 41, 84, 20, 85, 42, 86, 4, 87, 43, 88, 21, 89, 44, 90, 10, 91, 45,
  92, 22, 93, 46, 94, 0, 95, 47, 96, 23, 97, 48, 98, 11, 99,
  49, 100, 24, 101, 50, 102, 5, 103, 51, 104, 25, 105, 52, 106, 12, 107,
  53, 108, 26, 109, 54, 110, 2, 111, 55, 112, 27, 113, 56,
  114, 13, 115, 57, 116, 28, 117, 58, 118, 6, 119, 59, 120, 29, 121, 60,
  122, 14, 123, 61, 124, 30, 125, 62, 126, -2, 127, 63, 128,
  31, 129, 64, 130, 15, 131, 65, 132, 32, 133, 66, 134, 7, 135, 67, 136,
  33, 137, 68, 138, 16, 139, 69, 140, 34, 141, 70, 142, 3,
  143, 71, 144, 35, 145, 72, 146, 17, 147, 73, 148, 36, 149, 74, 150, 8,
  151, 75, 152, 37, 153, 76, 154, 18, 155, 77, 156, 38, 157,
  78, 158, 1, 159, 79, 160, 39, 161, 80, 162, 19, 163, 81, 164, 40, 165,
  82, 166, 9, 167, 83, 168, 41, 169, 84, 170, 20, 171, 85, 172,
  42, 173, 86, 174, 4, 175, 87, 176, 43, 177, 88, 178, 21, 179, 89, 180,
  44, 181, 90, 182, 10, 183, 91, 184, 45, 185, 92, 186, 22, 187,
  93, 188, 46, 189, 94, 190, 0, 191, 95, 192, 47, 193, 96, 194, 23, 195,
  97, 196, 48, 197, 98, 198, 11, 199, 99, 200, 49, 201, 100, 202,
  24, 203, 101, 204, 50, 205, 102, 206, 5, 207, 103, 208, 51, 209, 104,
  210, 25, 211, 105, 212, 52, 213, 106, 214, 12, 215, 107, 216, 53,
  217, 108, 218, 26, 219, 109, 220, 54, 221, 110, 222, 2, 223, 111, 224,
  55, 225, 112, 226, 27, 227, 113, 228, 56, 229, 114, 230, 13, 231,
  115, 232, 57, 233, 116, 234, 28, 235, 117, 236, 58, 237, 118, 238, 6,
  239, 119, 240, 59, 241, 120, 242, 29, 243, 121, 244, 60, 245, 122,
  246, 14, 247, 123, 248, 61, 249, 124, 250, 30, 251, 125, 252, 62, 253,
  126, 254, -2, 255, 127, 256, 63, 257, 128, 258, 31, 259, 129, 260,
  64, 261, 130, 262, 15, 263, 131, 264, 65, 265, 132, 266, 32, 267, 133,
  268, 66, 269, 134, 270, 7, 271, 135, 272, 67, 273, 136, 274, 33,
  275, 137, 276, 68, 277, 138, 278, 16, 279, 139, 280, 69, 281, 140, 282,
  34, 283, 141, 284, 70, 285, 142, 286, 3, 287, 143, 288, 71, 289,
  144, 290, 35, 291, 145, 292, 72, 293, 146, 294, 17, 295, 147, 296, 73,
  297, 148, 298, 36, 299, 149, 300, 74, 301, 150, 302, 8, 303, 151,
  304, 75, 305, 152, 306, 37, 307, 153, 308, 76, 309, 154, 310, 18, 311,
  155, 312, 77, 313, 156, 314, 38, 315, 157, 316, 78, 317, 158, 318,
  1, 319, 159, 320, 79, 321, 160, 322, 39, 323, 161, 324, 80, 325, 162,
  326, 19, 327, 163, 328, 81, 329, 164, 330, 40, 331, 165, 332, 82,
  333, 166, 334, 9, 335, 167, 336, 83, 337, 168, 338, 41, 339, 169, 340,
  84, 341, 170, 342, 20, 343, 171, 344, 85, 345, 172, 346, 42, 347,
  173, 348, 86, 349, 174, 350, 4, 351, 175, 352, 87, 353, 176, 354, 43,
  355, 177, 356, 88, 357, 178, 358, 21, 359, 179, 360, 89, 361, 180,
  362, 44, 363, 181, 364, 90, 365, 182, 366, 10, 367, 183, 368, 91, 369,
  184, 370, 45, 371, 185, 372, 92, 373, 186, 374, 22, 375, 187, 376,
  93, 377, 188, 378, 46, 379, 189, 380, 94, 381, 190, 382, 0, 383, 191,
  384, 95, 385, 192, 386, 47, 387, 193, 388, 96, 389, 194, 390, 23,
  391, 195, 392, 97, 393, 196, 394, 48, 395, 197, 396, 98, 397, 198, 398,
  11, 399, 199, 400, 99, 401, 200, 402, 49, 403, 201, 404, 100, 405,
  202, 406, 24, 407, 203, 408, 101, 409, 204, 410, 50, 411, 205, 412, 102,
  413, 206, 414, 5, 415, 207, 416, 103, 417, 208, 418, 51, 419,
  209, 420, 104, 421, 210, 422, 25, 423, 211, 424, 105, 425, 212, 426, 52,
  427, 213, 428, 106, 429, 214, 430, 12, 431, 215, 432, 107, 433,
  216, 434, 53, 435, 217, 436, 108, 437, 218, 438, 26, 439, 219, 440, 109,
  441, 220, 442, 54, 443, 221, 444, 110, 445, 222, 446, 2, 447,
  223, 448, 111, 449, 224, 450, 55, 451, 225, 452, 112, 453, 226, 454, 27,
  455, 227, 456, 113, 457, 228, 458, 56, 459, 229, 460, 114, 461,
  230, 462, 13, 463, 231, 464, 115, 465, 232, 466, 57, 467, 233, 468, 116,
  469, 234, 470, 28, 471, 235, 472, 117, 473, 236, 474, 58, 475,
  237, 476, 118, 477, 238, 478, 6, 479, 239, 480, 119, 481, 240, 482, 59,
  483, 241, 484, 120, 485, 242, 486, 29, 487, 243, 488, 121, 489,
  244, 490, 60, 491, 245, 492, 122, 493, 246, 494, 14, 495, 247, 496, 123,
  497, 248, 498, 61, 499, 249, 500, 124, 501, 250, 502, 30, 503,
  251, 504, 125, 505, 252, 506, 62, 507, 253, 508, 126, 509, 254, 510, -2,
  511, 255, 512, 127, 513, 256, 514, 63, 515, 257, 516, 128, 517,
  258, 518, 31, 519, 259, 520, 129, 521, 260, 522, 64, 523, 261, 524, 130,
  525, 262, 526, 15, 527, 263, 528, 131, 529, 264, 530, 65, 531,
  265, 532, 132, 533, 266, 534, 32, 535, 267, 536, 133, 537, 268, 538, 66,
  539, 269, 540, 134, 541, 270, 542, 7, 543, 271, 544, 135, 545,
  272, 546, 67, 547, 273, 548, 136, 549, 274, 550, 33, 551, 275, 552, 137,
  553, 276, 554, 68, 555, 277, 556, 138, 557, 278, 558, 16, 559,
  279, 560, 139, 561, 280, 562, 69, 563, 281, 564, 140, 565, 282, 566, 34,
  567, 283, 568, 141, 569, 284, 570, 70, 571, 285, 572, 142, 573,
  286, 574, 3, 575, 287, 576, 143, 577, 288, 578, 71, 579, 289, 580, 144,
  581, 290, 582, 35, 583, 291, 584, 145, 585, 292, 586, 72, 587,
  293, 588, 146, 589, 294, 590, 17, 591, 295, 592, 147, 593, 296, 594, 73,
  595, 297, 596, 148, 597, 298, 598, 36, 599, 299, 600, 149, 601,
  300, 602, 74, 603, 301, 604, 150, 605, 302, 606, 8, 607, 303, 608, 151,
  609, 304, 610, 75, 611, 305, 612, 152, 613, 306, 614, 37, 615,
  307, 616, 153, 617, 308, 618, 76, 619, 309, 620, 154, 621, 310, 622, 18,
  623, 311, 624, 155, 625, 312, 626, 77, 627, 313, 628, 156, 629,
  314, 630, 38, 631, 315, 632, 157, 633, 316, 634, 78, 635, 317, 636, 158,
  637, 318, 638, 1, 639, 319, 640, 159, 641, 320, 642, 79, 643,
  321, 644, 160, 645, 322, 646, 39, 647, 323, 648, 161, 649, 324, 650, 80,
  651, 325, 652, 162, 653, 326, 654, 19, 655, 327, 656, 163, 657,
  328, 658, 81, 659, 329, 660, 164, 661, 330, 662, 40, 663, 331, 664, 165,
  665, 332, 666, 82, 667, 333, 668, 166, 669, 334, 670, 9, 671,
  335, 672, 167, 673, 336, 674, 83, 675, 337, 676, 168, 677, 338, 678, 41,
  679, 339, 680, 169, 681, 340, 682, 84, 683, 341, 684, 170,
  685, 342, 686, 20, 687, 343, 688, 171, 689, 344, 690, 85, 691, 345, 692,
  172, 693, 346, 694, 42, 695, 347, 696, 173, 697, 348, 698,
  86, 699, 349, 700, 174, 701, 350, 702, 4, 703, 351, 704, 175, 705, 352,
  706, 87, 707, 353, 708, 176, 709, 354, 710, 43, 711, 355, 712,
  177, 713, 356, 714, 88, 715, 357, 716, 178, 717, 358, 718, 21, 719, 359,
  720, 179, 721, 360, 722, 89, 723, 361, 724, 180, 725, 362, 726,
  44, 727, 363, 728, 181, 729, 364, 730, 90, 731, 365, 732, 182, 733, 366,
  734, 10, 735, 367, 736, 183, 737, 368, 738, 91, 739, 369, 740,
  184, 741, 370, 742, 45, 743, 371, 744, 185, 745, 372, 746, 92, 747, 373,
  748, 186, 749, 374, 750, 22, 751, 375, 752, 187, 753, 376, 754,
  93, 755, 377, 756, 188, 757, 378, 758, 46, 759, 379, 760, 189, 761, 380,
  762, 94, 763, 381, 764, 190, 765, 382, 766, 0, 767, 383, 768, 191,
  769, 384, 770, 95, 771, 385, 772, 192, 773, 386, 774, 47, 775, 387, 776,
  193, 777, 388, 778, 96, 779, 389, 780, 194, 781, 390, 782, 23,
  783, 391, 784, 195, 785, 392, 786, 97, 787, 393, 788, 196, 789, 394,
  790, 48, 791, 395, 792, 197, 793, 396, 794, 98, 795, 397, 796, 198,
  797, 398, 798, 11, 799, 399, 800, 199, 801, 400, 802, 99, 803, 401,
  804, 200, 805, 402, 806, 49, 807, 403, 808, 201, 809, 404, 810, 100,
  811, 405, 812, 202, 813, 406, 814, 24, 815, 407, 816, 203, 817, 408,
  818, 101, 819, 409, 820, 204, 821, 410, 822, 50, 823, 411, 824, 205,
  825, 412, 826, 102, 827, 413, 828, 206, 829, 414, 830, 5, 831, 415,
  832, 207, 833, 416, 834, 103, 835, 417, 836, 208, 837, 418, 838, 51,
  839, 419, 840, 209, 841, 420, 842, 104, 843, 421, 844, 210, 845, 422,
  846, 25, 847, 423, 848, 211, 849, 424, 850, 105, 851, 425, 852, 212,
  853, 426, 854, 52, 855, 427, 856, 213, 857, 428, 858, 106, 859, 429,
  860, 214, 861, 430, 862, 12, 863, 431, 864, 215, 865, 432, 866, 107,
  867, 433, 868, 216, 869, 434, 870, 53, 871, 435, 872, 217, 873, 436,
  874, 108, 875, 437, 876, 218, 877, 438, 878, 26, 879, 439, 880, 219,
  881, 440, 882, 109, 883, 441, 884, 220, 885, 442, 886, 54, 887, 443,
  888, 221, 889, 444, 890, 110, 891, 445, 892, 222, 893, 446, 894, 2,
  895, 447, 896, 223, 897, 448, 898, 111, 899, 449, 900, 224, 901, 450,
  902, 55, 903, 451, 904, 225, 905, 452, 906, 112, 907, 453, 908, 226,
  909, 454, 910, 27, 911, 455, 912, 227, 913, 456, 914, 113, 915, 457,
  916, 228, 917, 458, 918, 56, 919, 459, 920, 229, 921, 460, 922, 114,
  923, 461, 924, 230, 925, 462, 926, 13, 927, 463, 928, 231, 929, 464,
  930, 115, 931, 465, 932, 232, 933, 466, 934, 57, 935, 467, 936, 233,
  937, 468, 938, 116, 939, 469, 940, 234, 941, 470, 942, 28, 943, 471,
  944, 235, 945, 472, 946, 117, 947, 473, 948, 236, 949, 474, 950, 58,
  951, 475, 952, 237, 953, 476, 954, 118, 955, 477, 956, 238, 957, 478,
  958, 6, 959, 479, 960, 239, 961, 480, 962, 119, 963, 481, 964, 240,
  965, 482, 966, 59, 967, 483, 968, 241, 969, 484, 970, 120, 971, 485,
  972, 242, 973, 486, 974, 29, 975, 487, 976, 243, 977, 488, 978, 121,
  979, 489, 980, 244, 981, 490, 982, 60, 983, 491, 984, 245, 985, 492,
  986, 122, 987, 493, 988, 246, 989, 494, 990, 14, 991, 495, 992, 247,
  993, 496, 994, 123, 995, 497, 996, 248, 997, 498, 998, 61, 999, 499,
  1000, 249, 1001, 500, 1002, 124, 1003, 501, 1004, 250, 1005, 502,
  1006, 30, 1007, 503, 1008, 251, 1009, 504, 1010, 125, 1011, 505, 1012,
  252, 1013, 506, 1014, 62, 1015, 507, 1016, 253, 1017, 508, 1018,
  126, 1019, 509, 1020, 254, 1021, 510, 1022
};


static int rightIdx[2047] = {
  -1, 0, -1, 1, 0, 2, -1, 3, 1, 4, 0, 5, 2, 6, -1, 7, 3, 8, 1, 9, 4, 10,
  0, 11, 5, 12, 2, 13, 6, 14, -1, 15, 7, 16,
  3, 17, 8, 18, 1, 19, 9, 20, 4, 21, 10, 22, 0, 23, 11, 24, 5, 25, 12,
  26, 2, 27, 13, 28, 6, 29, 14, 30, -1, 31, 15, 32, 7, 33, 16, 34, 3, 35,
  17, 36, 8, 37, 18, 38, 1, 39, 19, 40, 9, 41, 20, 42, 4, 43, 21, 44, 10,
  45, 22, 46, 0, 47, 23, 48, 11, 49, 24, 50, 5, 51, 25, 52, 12, 53, 26,
  54, 2, 55, 27, 56, 13, 57, 28, 58, 6, 59, 29, 60, 14, 61, 30, 62, -1,
  63, 31, 64, 15, 65, 32, 66, 7, 67, 33, 68, 16, 69, 34, 70, 3, 71, 35, 72,
  17, 73, 36, 74, 8, 75, 37, 76, 18, 77, 38, 78, 1, 79, 39, 80, 19, 81,
  40, 82, 9, 83, 41, 84, 20, 85, 42, 86, 4, 87, 43, 88, 21, 89, 44, 90, 10,
  91, 45, 92, 22, 93, 46, 94, 0, 95, 47, 96, 23, 97, 48, 98, 11, 99, 49,
  100, 24, 101, 50, 102, 5, 103, 51, 104, 25, 105, 52, 106, 12, 107, 53,
  108, 26, 109, 54, 110, 2, 111, 55, 112, 27, 113, 56, 114, 13, 115, 57,
  116, 28, 117, 58, 118, 6, 119, 59, 120, 29, 121, 60, 122, 14, 123, 61,
  124, 30, 125, 62, 126, -1, 127, 63, 128, 31, 129, 64, 130, 15, 131, 65,
  132, 32, 133, 66, 134, 7, 135, 67, 136, 33, 137, 68, 138, 16, 139, 69,
  140, 34, 141, 70, 142, 3, 143, 71, 144, 35, 145, 72, 146, 17, 147, 73,
  148, 36, 149, 74, 150, 8, 151, 75, 152, 37, 153, 76, 154, 18, 155, 77,
  156, 38, 157, 78, 158, 1, 159, 79, 160, 39, 161, 80, 162, 19, 163, 81,
  164, 40, 165, 82, 166, 9, 167, 83, 168, 41, 169, 84, 170, 20, 171, 85,
  172, 42, 173, 86, 174, 4, 175, 87, 176, 43, 177, 88, 178, 21, 179, 89,
  180, 44, 181, 90, 182, 10, 183, 91, 184, 45, 185, 92, 186, 22, 187, 93,
  188, 46, 189, 94, 190, 0, 191, 95, 192, 47, 193, 96, 194, 23, 195, 97,
  196, 48, 197, 98, 198, 11, 199, 99, 200, 49, 201, 100, 202, 24, 203, 101,
  204, 50, 205, 102, 206, 5, 207, 103, 208, 51, 209, 104, 210, 25, 211,
  105, 212, 52, 213, 106, 214, 12, 215, 107, 216, 53, 217, 108, 218, 26, 219,
  109, 220, 54, 221, 110, 222, 2, 223, 111, 224, 55, 225, 112, 226, 27,
  227, 113, 228, 56, 229, 114, 230, 13, 231, 115, 232, 57, 233, 116, 234, 28,
  235, 117, 236, 58, 237, 118, 238, 6, 239, 119, 240, 59, 241, 120, 242,
  29, 243, 121, 244, 60, 245, 122, 246, 14, 247, 123, 248, 61, 249, 124,
  250, 30, 251, 125, 252, 62, 253, 126, 254, -1, 255, 127, 256, 63, 257,
  128, 258, 31, 259, 129, 260, 64, 261, 130, 262, 15, 263, 131, 264, 65,
  265, 132, 266, 32, 267, 133, 268, 66, 269, 134, 270, 7, 271, 135, 272,
  67, 273, 136, 274, 33, 275, 137, 276, 68, 277, 138, 278, 16, 279, 139,
  280, 69, 281, 140, 282, 34, 283, 141, 284, 70, 285, 142, 286, 3, 287,
  143, 288, 71, 289, 144, 290, 35, 291, 145, 292, 72, 293, 146, 294, 17, 295,
  147, 296, 73, 297, 148, 298, 36, 299, 149, 300, 74, 301, 150, 302, 8,
  303, 151, 304, 75, 305, 152, 306, 37, 307, 153, 308, 76, 309, 154, 310, 18,
  311, 155, 312, 77, 313, 156, 314, 38, 315, 157, 316, 78, 317, 158, 318,
  1, 319, 159, 320, 79, 321, 160, 322, 39, 323, 161, 324, 80, 325, 162, 326,
  19, 327, 163, 328, 81, 329, 164, 330, 40, 331, 165, 332, 82, 333, 166,
  334, 9, 335, 167, 336, 83, 337, 168, 338, 41, 339, 169, 340, 84, 341, 170,
  342, 20, 343, 171, 344, 85, 345, 172, 346, 42, 347, 173, 348, 86, 349,
  174, 350, 4, 351, 175, 352, 87, 353, 176, 354, 43, 355, 177, 356, 88, 357,
  178, 358, 21, 359, 179, 360, 89, 361, 180, 362, 44, 363, 181, 364, 90,
  365, 182, 366, 10, 367, 183, 368, 91, 369, 184, 370, 45, 371, 185, 372, 92,
  373, 186, 374, 22, 375, 187, 376, 93, 377, 188, 378, 46, 379, 189, 380,
  94, 381, 190, 382, 0, 383, 191, 384, 95, 385, 192, 386, 47, 387, 193, 388,
  96, 389, 194, 390, 23, 391, 195, 392, 97, 393, 196, 394, 48, 395, 197,
  396, 98, 397, 198, 398, 11, 399, 199, 400, 99, 401, 200, 402, 49, 403, 201,
  404, 100, 405, 202, 406, 24, 407, 203, 408, 101, 409, 204, 410, 50, 411,
  205, 412, 102, 413, 206, 414, 5, 415, 207, 416, 103, 417, 208, 418, 51,
  419, 209, 420, 104, 421, 210, 422, 25, 423, 211, 424, 105, 425, 212, 426,
  52, 427, 213, 428, 106, 429, 214, 430, 12, 431, 215, 432, 107, 433, 216,
  434, 53, 435, 217, 436, 108, 437, 218, 438, 26, 439, 219, 440, 109, 441,
  220, 442, 54, 443, 221, 444, 110, 445, 222, 446, 2, 447, 223, 448, 111,
  449, 224, 450, 55, 451, 225, 452, 112, 453, 226, 454, 27, 455, 227, 456,
  113, 457, 228, 458, 56, 459, 229, 460, 114, 461, 230, 462, 13, 463, 231,
  464, 115, 465, 232, 466, 57, 467, 233, 468, 116, 469, 234, 470, 28, 471,
  235, 472, 117, 473, 236, 474, 58, 475, 237, 476, 118, 477, 238, 478, 6,
  479, 239, 480, 119, 481, 240, 482, 59, 483, 241, 484, 120, 485, 242, 486,
  29, 487, 243, 488, 121, 489, 244, 490, 60, 491, 245, 492, 122, 493, 246,
  494, 14, 495, 247, 496, 123, 497, 248, 498, 61, 499, 249, 500, 124, 501,
  250, 502, 30, 503, 251, 504, 125, 505, 252, 506, 62, 507, 253, 508, 126,
  509, 254, 510, -1, 511, 255, 512, 127, 513, 256, 514, 63, 515, 257, 516,
  128, 517, 258, 518, 31, 519, 259, 520, 129, 521, 260, 522, 64, 523, 261,
  524, 130, 525, 262, 526, 15, 527, 263, 528, 131, 529, 264, 530, 65, 531,
  265, 532, 132, 533, 266, 534, 32, 535, 267, 536, 133, 537, 268, 538, 66,
  539, 269, 540, 134, 541, 270, 542, 7, 543, 271, 544, 135, 545, 272, 546,
  67, 547, 273, 548, 136, 549, 274, 550, 33, 551, 275, 552, 137, 553, 276,
  554, 68, 555, 277, 556, 138, 557, 278, 558, 16, 559, 279, 560, 139, 561,
  280, 562, 69, 563, 281, 564, 140, 565, 282, 566, 34, 567, 283, 568, 141,
  569, 284, 570, 70, 571, 285, 572, 142, 573, 286, 574, 3, 575, 287, 576,
  143, 577, 288, 578, 71, 579, 289, 580, 144, 581, 290, 582, 35, 583, 291,
  584, 145, 585, 292, 586, 72, 587, 293, 588, 146, 589, 294, 590, 17, 591,
  295, 592, 147, 593, 296, 594, 73, 595, 297, 596, 148, 597, 298, 598, 36,
  599, 299, 600, 149, 601, 300, 602, 74, 603, 301, 604, 150, 605, 302, 606,
  8, 607, 303, 608, 151, 609, 304, 610, 75, 611, 305, 612, 152, 613, 306,
  614, 37, 615, 307, 616, 153, 617, 308, 618, 76, 619, 309, 620, 154, 621,
  310, 622, 18, 623, 311, 624, 155, 625, 312, 626, 77, 627, 313, 628, 156,
  629, 314, 630, 38, 631, 315, 632, 157, 633, 316, 634, 78, 635, 317, 636,
  158, 637, 318, 638, 1, 639, 319, 640, 159, 641, 320, 642, 79, 643, 321,
  644, 160, 645, 322, 646, 39, 647, 323, 648, 161, 649, 324, 650, 80, 651,
  325, 652, 162, 653, 326, 654, 19, 655, 327, 656, 163, 657, 328, 658, 81,
  659, 329, 660, 164, 661, 330, 662, 40, 663, 331, 664, 165, 665, 332, 666,
  82, 667, 333, 668, 166, 669, 334, 670, 9, 671, 335, 672, 167, 673, 336,
  674, 83, 675, 337, 676, 168, 677, 338, 678, 41, 679, 339, 680, 169, 681,
  340, 682, 84, 683, 341, 684, 170, 685, 342, 686, 20, 687, 343, 688, 171,
  689, 344, 690, 85, 691, 345, 692, 172, 693, 346, 694, 42, 695, 347, 696,
  173, 697, 348, 698, 86, 699, 349, 700, 174, 701, 350, 702, 4, 703, 351,
  704, 175, 705, 352, 706, 87, 707, 353, 708, 176, 709, 354, 710, 43, 711,
  355, 712, 177, 713, 356, 714, 88, 715, 357, 716, 178, 717, 358, 718, 21,
  719, 359, 720, 179, 721, 360, 722, 89, 723, 361, 724, 180, 725, 362, 726,
  44, 727, 363, 728, 181, 729, 364, 730, 90, 731, 365, 732, 182, 733, 366,
  734, 10, 735, 367, 736, 183, 737, 368, 738, 91, 739, 369, 740, 184, 741,
  370, 742, 45, 743, 371, 744, 185, 745, 372, 746, 92, 747, 373, 748, 186,
  749, 374, 750, 22, 751, 375, 752, 187, 753, 376, 754, 93, 755, 377, 756,
  188, 757, 378, 758, 46, 759, 379, 760, 189, 761, 380, 762, 94, 763, 381,
  764, 190, 765, 382, 766, 0, 767, 383, 768, 191, 769, 384, 770, 95, 771,
  385, 772, 192, 773, 386, 774, 47, 775, 387, 776, 193, 777, 388, 778, 96,
  779, 389, 780, 194, 781, 390, 782, 23, 783, 391, 784, 195, 785, 392, 786,
  97, 787, 393, 788, 196, 789, 394, 790, 48, 791, 395, 792, 197, 793, 396,
  794, 98, 795, 397, 796, 198, 797, 398, 798, 11, 799, 399, 800, 199, 801,
  400, 802, 99, 803, 401, 804, 200, 805, 402, 806, 49, 807, 403, 808, 201,
  809, 404, 810, 100, 811, 405, 812, 202, 813, 406, 814, 24, 815, 407, 816,
  203, 817, 408, 818, 101, 819, 409, 820, 204, 821, 410, 822, 50, 823, 411,
  824, 205, 825, 412, 826, 102, 827, 413, 828, 206, 829, 414, 830, 5, 831,
  415, 832, 207, 833, 416, 834, 103, 835, 417, 836, 208, 837, 418, 838, 51,
  839, 419, 840, 209, 841, 420, 842, 104, 843, 421, 844, 210, 845, 422, 846,
  25, 847, 423, 848, 211, 849, 424, 850, 105, 851, 425, 852, 212, 853, 426,
  854, 52, 855, 427, 856, 213, 857, 428, 858, 106, 859, 429, 860, 214, 861,
  430, 862, 12, 863, 431, 864, 215, 865, 432, 866, 107, 867, 433, 868, 216,
  869, 434, 870, 53, 871, 435, 872, 217, 873, 436, 874, 108, 875, 437, 876,
  218, 877, 438, 878, 26, 879, 439, 880, 219, 881, 440, 882, 109, 883, 441,
  884, 220, 885, 442, 886, 54, 887, 443, 888, 221, 889, 444, 890, 110, 891,
  445, 892, 222, 893, 446, 894, 2, 895, 447, 896, 223, 897, 448, 898, 111,
  899, 449, 900, 224, 901, 450, 902, 55, 903, 451, 904, 225, 905, 452, 906,
  112, 907, 453, 908, 226, 909, 454, 910, 27, 911, 455, 912, 227, 913, 456,
  914, 113, 915, 457, 916, 228, 917, 458, 918, 56, 919, 459, 920, 229, 921,
  460, 922, 114, 923, 461, 924, 230, 925, 462, 926, 13, 927, 463, 928, 231,
  929, 464, 930, 115, 931, 465, 932, 232, 933, 466, 934, 57, 935, 467, 936,
  233, 937, 468, 938, 116, 939, 469, 940, 234, 941, 470, 942, 28, 943, 471,
  944, 235, 945, 472, 946, 117, 947, 473, 948, 236, 949, 474, 950, 58, 951,
  475, 952, 237, 953, 476, 954, 118, 955, 477, 956, 238, 957, 478, 958, 6,
  959, 479, 960, 239, 961, 480, 962, 119, 963, 481, 964, 240, 965, 482, 966,
  59, 967, 483, 968, 241, 969, 484, 970, 120, 971, 485, 972, 242, 973, 486,
  974, 29, 975, 487, 976, 243, 977, 488, 978, 121, 979, 489, 980, 244, 981,
  490, 982, 60, 983, 491, 984, 245, 985, 492, 986, 122, 987, 493, 988, 246,
  989, 494, 990, 14, 991, 495, 992, 247, 993, 496, 994, 123, 995, 497, 996,
  248, 997, 498, 998, 61, 999, 499, 1000, 249, 1001, 500, 1002, 124, 1003,
  501, 1004, 250, 1005, 502, 1006, 30, 1007, 503, 1008, 251, 1009, 504, 1010,
  125, 1011, 505, 1012, 252, 1013, 506, 1014, 62, 1015, 507, 1016, 253,
  1017, 508, 1018, 126, 1019, 509, 1020, 254, 1021, 510, 1022, -1
};

Stretching::Stretching(size_t dimension, const BoundingBox1D* boundaries,
                      const Stretching1D* t) : BoundingBox(dimension, boundaries) {
  bTrivialCube = true;
  stretching1Ds = new Stretching1D[dimension];
  dimensionBoundaries = new BoundingBox1D[dimension];
  discreteVectorLevel = new int[dimension];
  stretchingMode = new std::string("analytic");

  for (size_t i = 0; i < dimension; i++) {
    dimensionBoundaries[i] = boundaries[i];

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }

    discreteVectorLevel[i] = -1;
    stretching1Ds[i] = t[i];
  }

  generateLookupTable();
}

Stretching::Stretching(size_t dimension, const std::vector<BoundingBox1D>& boundaries,
                       const std::vector<Stretching1D>& t) : BoundingBox(dimension) {
  bTrivialCube = true;
  stretching1Ds = new Stretching1D[dimension];
  dimensionBoundaries = new BoundingBox1D[dimension];
  discreteVectorLevel = new int[dimension];
  stretchingMode = new std::string("analytic");

  for (size_t i = 0; i < dimension; i++) {
    dimensionBoundaries[i] = boundaries[i];

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }

    discreteVectorLevel[i] = -1;
    stretching1Ds[i] = t[i];
  }

  generateLookupTable();
}

Stretching::Stretching(size_t dimension, std::vector<double>* coordinates) :
    BoundingBox(dimension) {
  bTrivialCube = true;
  stretching1Ds = new Stretching1D[dimension];
  dimensionBoundaries = new BoundingBox1D[dimension];
  discreteVectorLevel = new int[dimension];
  stretchingMode = new std::string("discrete");
  int _1DArrayLength = 0;

  for (size_t i = 0; i < dimension; i++) {
    _1DArrayLength = static_cast<int>(coordinates[i].size());
    dimensionBoundaries[i].leftBoundary = coordinates[i][0];
    dimensionBoundaries[i].rightBoundary = coordinates[i][_1DArrayLength - 1];
    dimensionBoundaries[i].bDirichletLeft = true;
    dimensionBoundaries[i].bDirichletRight = true;

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }

    parseVectorToLookupTable(coordinates[i], stretching1Ds[i], i,
                             discreteVectorLevel[i]);
    generateLeftRightArrays(stretching1Ds[i], i);
  }

  // No need to generate lookup table,
  // as it is done in parseVectorToLookupTable
}

Stretching::Stretching(const Stretching& copyStretching) : BoundingBox(copyStretching) {
  bTrivialCube = true;
  dimensionBoundaries = new BoundingBox1D[dimension];
  stretching1Ds = new Stretching1D[dimension];
  discreteVectorLevel = new int[dimension];
  stretchingMode = new std::string(*(copyStretching.getStretchingMode()));


  for (size_t i = 0; i < dimension; i++) {
    dimensionBoundaries[i] = copyStretching.getBoundary(i);

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }

    discreteVectorLevel[i] = copyStretching.discreteVectorLevel[i];
    stretching1Ds[i] =
      copyStretching.getStretching1D(*reinterpret_cast<int*>(&i));
  }
}

Stretching::~Stretching() {
  if (stretching1Ds != NULL) {
    delete[] stretching1Ds;
  }

  if (discreteVectorLevel != NULL) {
    delete[] discreteVectorLevel;
  }

  if (stretchingMode != NULL) {
    delete stretchingMode;
  }
}


void Stretching::generateLookupTable() {
  //  std::cout<<"Generate Lookup Table"<<std::endl;
  for (size_t i = 0; i < dimension; i++) {
    if (stretching1Ds[i].type == "log") {
      logXform(stretching1Ds[i], i);
    } else if (stretching1Ds[i].type == "sinh") {
      leentvaarXform(stretching1Ds[i], i);
    } else if (stretching1Ds[i].type == "fitob") {
      /*
       * This means we copied this type using a copy constructor,
       *  no need for calculation necessary
       *  Should not occur in a normal operation.
       */
    } else {
      // std::cout << "Stretching::generateLookupTable :
      // analytic stretching type not supported!" << std::endl;
      // throw application_exception("Stretching::generateLookupTable :
      // analytic stretching type not supported!");
      noXform(stretching1Ds[i], i);  // in case of discrete stretching
    }

    generateLeftRightArrays(stretching1Ds[i], i);
  }
}

void Stretching::generateLeftRightArrays(Stretching1D& str1D, size_t d) {
  int idx = 0;
  int lidx = 0, ridx = 0;

  for (int l = 1; l <= (LOOKUPMAX); l++) {
    int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

    for (int i = 1; i <= elemPerLevel; i++) {
      // Create left neighbors table
      lidx = leftIdx[idx];

      if (lidx >= 0) {
        str1D.lookup[idx][1] = str1D.lookup[lidx][0];
      } else if (lidx == -2) {
        str1D.lookup[idx][1] =  dimensionBoundaries[d].leftBoundary;
      }

      // Create right neighbors table
      ridx = rightIdx[idx];

      if (ridx >= 0) {
        str1D.lookup[idx][2] = str1D.lookup[ridx][0];
      } else if (ridx == -1) {
        str1D.lookup[idx][2] =  dimensionBoundaries[d].rightBoundary;
      }

      idx++;
    }
  }
}

double Stretching::stretchingXform(int level, int index, size_t d) const {
  double temp = -1;
  // must be updated!!!!!!

  std::string t = stretching1Ds[d].type;

  // refer to the lookup table
  if (level <= LOOKUPMAX) {
    temp = stretching1Ds[d].lookup[calculateLookupIndex(level, index)][0];
  } else if (t == "log") {
    // calculate using the appropriate function
    temp = logXform(level, index, d);
  } else if (t == "sinh") {
    temp = leentvaarXform(level, index, d);
  } else if (t == "fitob") {
    int pow2deltaL = static_cast<int>(pow(2.0, (level - LOOKUPMAX)));
    double dIndex = static_cast<double>(index) /
                     static_cast<double>(pow2deltaL);
    bool leftContinue = false, rightContinue = false;
    int rightLevel = LOOKUPMAX, leftLevel = LOOKUPMAX;
    int leftIndex = static_cast<int>(floor(dIndex));
    int rightIndex = static_cast<int>(ceil(dIndex));

    if (leftIndex % 2 == 0) {
      leftContinue = true;
    }

    if (rightIndex % 2 == 0) {
      rightContinue = true;
    }

    bool loopContinue = leftContinue || rightContinue;

    while (loopContinue) {
      if (leftContinue) {
        leftIndex = leftIndex / 2;
        leftLevel = leftLevel - 1;

        // if left index odd, or left level 0, stop calculating left part
        if (((leftIndex % 2 != 0)) || (leftLevel == 0)) {
          if (leftLevel == 0) {
            leftIndex = 0;
          }

          leftContinue = false;
        }
      }

      if (rightContinue) {
        rightIndex = rightIndex / 2;
        rightLevel = rightLevel - 1;

        // if right index odd, or right level 0, stop calculating right part
        if (((rightIndex % 2) != 0) || (rightLevel == 0)) {
          if (rightLevel == 0) {
            rightIndex = 1;
          }

          rightContinue = false;
        }
      }

      loopContinue = rightContinue || leftContinue;
    }

    double posl = getCoordinate(leftLevel, leftIndex, d);
    double posr = getCoordinate(rightLevel, rightIndex, d);
    double step = (posr - posl) / static_cast<double>(pow2deltaL);
    temp = posl + step * (dIndex - floor(dIndex)) * pow2deltaL;
  } else {
    temp = noXform(level, index, d);
  }

  return temp;
}

void Stretching::logXform(Stretching1D& str1D, size_t d) const {
  int idx = 0;
  double f_a = log(getBoundary(d).leftBoundary);
  double f_b = log(getBoundary(d).rightBoundary);

  for (int l = 1; l <= (LOOKUPMAX); l++) {
    int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

    for (int i = 1; i <= elemPerLevel; i++) {
      str1D.lookup[idx++][0] =
        exp(f_a + static_cast<double>(2 * i - 1) /
            static_cast<double>(2 * elemPerLevel) * (f_b - f_a));
    }
  }
}

double Stretching::logXform(int l, int i, size_t d) const {
  double f_a = log(getBoundary(d).leftBoundary);
  double f_b = log(getBoundary(d).rightBoundary);
  int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

  return exp (f_a + static_cast<double>(i) /
              static_cast<double>(2 * elemPerLevel) * (f_b - f_a));
}

void Stretching::leentvaarXform(Stretching1D& str1D, size_t d) const {
  // f_inv(sinh(x))= log(x+sqrt(x²+1))
  int idx = 0;
  double a = (getBoundary(d).leftBoundary);
  double b = (getBoundary(d).rightBoundary);
  double sinhArgumenta = ((a - str1D.x_0) * str1D.xsi);
  double sinhArgumentb = ((b - str1D.x_0) * str1D.xsi);
  double f_a = log(sinhArgumenta + sqrt(sinhArgumenta * sinhArgumenta + 1));
  double f_b = log(sinhArgumentb + sqrt(sinhArgumentb * sinhArgumentb + 1));

  for (int l = 1; l <= (LOOKUPMAX); l++) {
    int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

    for (int i = 1; i <= elemPerLevel; i++) {
      str1D.lookup[idx++][0] =
        (1 / str1D.xsi * sinh((f_a + static_cast<double>(2 * i - 1) /
                               static_cast<double>(2 * elemPerLevel) *
                               (f_b - f_a))) + str1D.x_0);
    }
  }
}

double Stretching::leentvaarXform(int l, int i, size_t d) const {
  double a = (getBoundary(d).leftBoundary);
  double b = (getBoundary(d).rightBoundary);
  double sinhArgumenta = ((a - stretching1Ds[d].x_0) *
                           stretching1Ds[d].xsi);
  double sinhArgumentb = ((b - stretching1Ds[d].x_0) *
                           stretching1Ds[d].xsi);
  double f_a = log(sinhArgumenta + sqrt(sinhArgumenta * sinhArgumenta + 1));
  double f_b = log(sinhArgumentb + sqrt(sinhArgumentb * sinhArgumentb + 1));
  int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

  return (1 / stretching1Ds[d].xsi * sinh((f_a + static_cast<double>(i) /
          static_cast<double>(2 * elemPerLevel) *
          (f_b - f_a))) + stretching1Ds[d].x_0);
}

void Stretching::noXform(Stretching1D& str1D, size_t d) const {
  int idx = 0;
  double f_a = getBoundary(d).leftBoundary;
  double f_b = getBoundary(d).rightBoundary;

  for (int l = 1; l <= (LOOKUPMAX); l++) {
    int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

    for (int i = 1; i <= elemPerLevel; i++) {
      str1D.lookup[idx++][0] =
        (f_a + static_cast<double>(2 * i - 1) /
         static_cast<double>(2 * elemPerLevel) * (f_b - f_a));
    }
  }
}

double Stretching::noXform(int l, int i, size_t dim) const {
  double f_a = getBoundary(dim).leftBoundary;
  double f_b = getBoundary(dim).rightBoundary;
  int elemPerLevel = static_cast<int>(pow(2.0, l - 1));

  return (f_a + static_cast<double>(i) /
          static_cast<double>(2 * elemPerLevel) * (f_b - f_a));
}



int Stretching::calculateLookupIndex(int l, int i) const {
  // Note: indices are always odd ;)
  return static_cast<int>((pow(2.0, l - 1) - 1 + (i - 1) / 2));
}

double Stretching::getCoordinate(int level, int index, size_t d) const {
  if (level == 0) {
    if (index == 0) {
      return dimensionBoundaries[d].leftBoundary;
    } else if (index == 1) {
      return dimensionBoundaries[d].rightBoundary;
    }
  } else {
    return stretchingXform(level, index, d);
  }

  // should not happen
  return 0.0;
}

Stretching1D Stretching::getStretching1D(size_t d) const {
  return stretching1Ds[d];
}

void Stretching::printLookupTable() const {
  for (size_t i = 0; i < dimension; i++) {
    std::cout << std::endl << "dim" << i << std::endl;

    for (int j = 0; j < LOOKUPSIZE; j++) {
      std::cout << stretching1Ds[i].lookup[j][0] << " ";
    }

    std::cout << std::endl;
  }
}

void Stretching::getAdjacentPositions(int level, int index, size_t d,
                                      double& posc, double& posl,
                                      double& posr) const {
  if (level <= LOOKUPMAX) {
    int idx = calculateLookupIndex(level, index);
    posc = stretching1Ds[d].lookup[idx][0];
    posl = stretching1Ds[d].lookup[idx][1];
    posr = stretching1Ds[d].lookup[idx][2];
  } else {
    int leftIndex, rightIndex, rightLevel, leftLevel;
    calculateNeighborSpecs(level, index, leftLevel, leftIndex, rightLevel,
                           rightIndex);
    posl = getCoordinate(leftLevel, leftIndex, d);
    posr = getCoordinate(rightLevel, rightIndex, d);
    posc = getCoordinate(level, index, d);
  }
}

void Stretching::calculateNeighborSpecs(int level, int index, int& leftLevel,
                                        int& leftIndex, int& rightLevel,
                                        int& rightIndex) const {
  bool loopContinue = true;
  bool leftContinue = true, rightContinue = true;

  leftIndex = index - 1;
  rightIndex = index + 1;
  rightLevel = level;
  leftLevel = level;

  if (leftIndex == 0) {
    leftLevel = 0;
    leftContinue = false;
  }

  if (rightIndex == pow(2.0, level)) {
    rightLevel = 0;
    rightIndex = 1;
    rightContinue = false;
  }

  loopContinue = leftContinue || rightContinue;

  while (loopContinue) {
    if (leftContinue) {
      leftIndex = leftIndex / 2;
      leftLevel = leftLevel - 1;

      // if left index odd, or left level 0, stop calculating left part
      if (((leftIndex % 2 != 0)) || (leftLevel == 0)) {
        if (leftLevel == 0) {
          leftIndex = 0;
        }

        leftContinue = false;
      }
    }

    if (rightContinue) {
      rightIndex = rightIndex / 2;
      rightLevel = rightLevel - 1;

      // if right index odd, or right level 0, stop calculating right part
      if (((rightIndex % 2) == 1) || (rightLevel == 0)) {
        if (rightLevel == 0) {
          rightIndex = 1;
        }

        rightContinue = false;
      }
    }

    loopContinue = rightContinue || leftContinue;
  }
}
void Stretching::parseVectorToLookupTable(std::vector<double>& vec,
                                          Stretching1D& stretch1d,
                                          size_t d, int& discreteVectorLevel) const {
  int level = static_cast<int>(log(static_cast<double>(vec.size() - 1)) / log(
                                 2.0));  // log2(vec.size()-1);

  if ((static_cast<size_t>(1) << level) + 1 != vec.size()) {
    std::cout << "parseVectorToLookupTable: "
              "Vector Size does not match, should be 2^l+1\n";
    return;
  }

  int rLevel, rIndex, lLevel, lIndex;
  discreteVectorLevel = level;
  int elemPerLevel = 0, idx = 0;
  double posl = 0, posr = 0;

  stretch1d.type = "fitob";

  // Handle current and lower levels
  for (int l = level; l > 0; l--) {
    elemPerLevel = static_cast<int>(pow(2.0, l));

    for (int i = 1; i < elemPerLevel; i = i + 2) {
      stretch1d.lookup[calculateLookupIndex(l, i)][0] = vec[i];
      vec[idx++] = vec[i - 1];
    }

    vec[idx++] = vec[elemPerLevel];
    idx = 0;
  }

  // Handle upper levels (up to LOOKUPMAX)
  // Note: Takes a long time due to getAdjacentPositions
  for (int l = level + 1; l <= LOOKUPMAX; l++) {
    elemPerLevel = static_cast<int>(pow(2.0, l));

    for (int i = 1; i < elemPerLevel; i = i + 2) {
      calculateNeighborSpecs(l, i, lLevel, lIndex, rLevel, rIndex);
      posl = getCoordinate(lLevel, lIndex, d);
      posr = getCoordinate(rLevel, rIndex, d);
      stretch1d.lookup[calculateLookupIndex(l, i)][0] =
        posl + 0.5 * (posr - posl);
    }
  }
}

std::string* Stretching::getStretchingMode() const {
  return stretchingMode;
}

std::vector<double>* Stretching::getDiscreteVector(bool bSort) const {
  std::vector<double>* vec;
  int elemsToRead = 0;

  if (*stretchingMode == "discrete") {
    vec = new std::vector<double>[dimension];

    for (size_t i = 0; i < dimension; i++) {
      Stretching1D str1d = getStretching1D(*reinterpret_cast<int*>(&i));
      elemsToRead = static_cast<int>((pow(2.0, discreteVectorLevel[i]) - 1));
      vec[i] = std::vector<double>(elemsToRead + 2, 0);
      vec[i][0] = dimensionBoundaries[i].leftBoundary;

      for (int j = 0; j < elemsToRead; j++) {
        vec[i][j + 1] = str1d.lookup[j][0];
      }

      vec[i][elemsToRead + 1] = dimensionBoundaries[i].rightBoundary;
    }

    if (bSort) {
      for (size_t i = 0; i < dimension; i++) {
        std::sort(vec[i].begin(), vec[i].end());
      }
    }

    return vec;
  } else {
    std::cout << "Unsupported stretchingMode\n";
    return nullptr;
  }
}

int* Stretching::getDiscreteVectorLevel() const {
  return discreteVectorLevel;
}

void Stretching::calculateNeighborLookup(int maxlevel) const {
  //  std::string file1 = "leftIndex.txt";
  //  std::string file2 = "rightIndex.txt";
  std::ofstream outfile1;
  std::ofstream outfile2;

  outfile1.open("leftIndex.txt");
  outfile2.open("rightIndex.txt");

  outfile1 << "{";
  outfile2 << "{";

  int lLevel, lIndex, rLevel, rIndex, temp;
  temp = 0;
  double elemPerLevel;

  for (int l = 1; l <= maxlevel; l++) {
    elemPerLevel = pow(2.0, l);

    for (int i = 1; i < elemPerLevel; i = i + 2) {
      calculateNeighborSpecs(l, i, lLevel, lIndex, rLevel, rIndex);

      // Left position
      if (lLevel == 0) {
        if (lIndex == 1)
          temp = -1;

        if (lIndex == 0)
          temp = -2;
      } else {
        temp = static_cast<int>((pow(2.0, lLevel - 1) - 1 + (lIndex - 1) / 2));
      }

      outfile1 << temp;
      outfile1 << ", ";

      // Right position
      if (rLevel == 0) {
        if (rIndex == 1)
          temp = -1;

        if (rIndex == 0)
          temp = -2;
      } else {
        temp = static_cast<int>((pow(2.0, rLevel - 1) - 1 + (rIndex - 1) / 2));
      }

      outfile2 << temp;
      outfile2 << ", ";
    }
  }

  outfile1 << "}";
  outfile2 << "}";

  outfile1.close();
  outfile2.close();
}

}  // namespace base
}  // namespace sgpp
