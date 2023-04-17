import numpy as np
import pandas as pd
import sys

def get_the_float(feature_value, number=8):
    float_list = [float(x) for x in feature_value]
    float_array = np.array(float_list)
    return list(np.round(float_array, number))

class PCP9:
    def __init__(self):
        pass

    def get_feature_name(self, se_len):
        p_property = ['Hydrophobicity', 'Hydrophilicity', 'Mass', 'pK1', 'pK2', 'pI', 'Rigidity', 'Flexibility', 'Irreplaceability']
        feature_name_list = []
        for el in p_property:
            for i in range(1, se_len + 1):
                feature_name_temp = 'PCP9_' + el + '_pos_' + str(i)
                feature_name_list.append(feature_name_temp)
        return feature_name_list

    def main(self, sequence):
        f_PCP = pd.read_csv('nine_physicochemical_properties_of_amino_acid_Stand.txt', sep='\t')
        AA_char = 'ACDEFGHIKLMNPQRSTVWY'
        feature_temp_dict = [dict() for i in range(9)]
        for char in AA_char:
            for i in range(9):
                feature_temp_dict[i][char] = f_PCP.loc[i][char]
        feature_temp_list = [list() for i in range(9)]
        for char in sequence:
            if char != 'X':
                for i in range(9):
                    feature_temp_list[i].append(str(feature_temp_dict[i][char]))
            else:
                for i in range(9):
                    feature_temp_list[i].append('0')
        feature_total = []
        for i in range(9):
            feature_total.extend(feature_temp_list[i])

        feature = get_the_float(feature_total)
        return feature

class PWAA:
    def __init__(self):
        self.AA_list = 'ACDEFGHIKLMNPQRSTVWY'

    def get_feature_name(self):
        feature_name_list = []
        for element in self.AA_list:
            featureName = 'PWAA_feature_%s' % element
            feature_name_list.append(featureName)
        return feature_name_list

    def main(self, sequence):
        length_up_down = (len(sequence) - 1) / 2
        feature = []
        for aa_char in self.AA_list:
            sum_inter = 0
            if aa_char not in sequence:
                feature.append(0)
            else:
                for sequence_index, sequence_char in enumerate(sequence):
                    if sequence_char == aa_char:
                        j = sequence_index - length_up_down
                        sum_inter = sum_inter + (j + abs(j) / length_up_down)
                c = (1 / (length_up_down * (length_up_down + 1))) * sum_inter
                feature.append(c)
        feature_total = get_the_float(feature)
        return feature_total

class ARPC:
    def __init__(self):
        self.AA_dict = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
                        'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, "X": 0}

    def get_feature_name(self, sequence_length):
        feature_name_list = []
        up_down_length = (sequence_length - 1) / 2
        for i in range(0, sequence_length):
            temp_site = i - up_down_length
            featureName = 'ARPC_site_%d' % int(temp_site)
            feature_name_list.append(featureName)
        return feature_name_list

    def main(self, sequence):
        length_up_down = (len(sequence) - 1) / 2
        feature = []
        for i in range(0, len(sequence)):
            site_value = (i - length_up_down)
            acid_code = self.AA_dict[sequence[i]]
            temp_feature = site_value * acid_code
            feature.append(temp_feature)
        feature_total = get_the_float(feature)
        return feature_total

class K_space:
    def __init__(self, k):
        self.k = k
        self.feature_name = []
        self.AA_list_sort = 'ACDEFGHIKLMNPQRSTVWY'

    def get_feature_name(self):
        for el in range(self.k + 1):
            prefix = 'Kspace_%d_space_' % el
            for i1 in self.AA_list_sort:
                for j1 in self.AA_list_sort:
                    self.feature_name.append(prefix + i1 + j1)
        return self.feature_name

    def main(self, sequence):
        two_mer_name = []
        for i2 in self.AA_list_sort:
            for j2 in self.AA_list_sort:
                combination = i2 + j2
                two_mer_name.append(combination)
        two_mer_dict = {}
        for item in two_mer_name:
            two_mer_dict[item] = 0
        for i, v in enumerate(sequence):
            if i + 2 + self.k > len(sequence):
                break
            else:
                now_k_mer = sequence[i:i + 2 + self.k:self.k + 1]
                if now_k_mer in two_mer_dict:
                    two_mer_dict[now_k_mer] += 1
        feature_vector = []
        for key, value in two_mer_dict.items():
            frequency = value / (len(sequence) - 1 - self.k)
            feature_vector.append(frequency)
        feature = get_the_float(feature_vector)
        return feature

class EBGW:
    def __init__(self, number_sub_sequence):
        self.number_sub_sequence = number_sub_sequence

    def get_feature_name(self):
        feature_name = []
        for i in range(1, self.number_sub_sequence + 1):
            prefix = 'EBGW_%d_feature_' % i
            total_number_sub_sequence = 3 * i
            for j in range(1, total_number_sub_sequence + 1):
                featureName = prefix + str(j)
                feature_name.append(featureName)
        return feature_name

    def main(self, sequence):
        cs = ['AFGILMPVW','CNQSTY','KHR','DE']
        hs_bi_sequence = ['','','']

        for i in range(3):
            for element in sequence:
                if element in cs[0] or element in cs[i+1]:
                    hs_bi_sequence[i] += str(1)
                else:
                    hs_bi_sequence[i] += str(0)
        feature = []
        feature_Hs_list = [list() for i in range(3)]

        for i in range(1, self.number_sub_sequence + 1):
            i_length = round(i * len(sequence) / self.number_sub_sequence)
            for j in range(3):
                feature_hs_sub_sequence = hs_bi_sequence[j][0:i_length].count('1') / i_length
                feature_Hs_list[j].append(feature_hs_sub_sequence)

        for feature_Hs in feature_Hs_list:
            for el in feature_Hs:
                feature.append(el)

        feature_total = get_the_float(feature)
        return feature_total

class CTD:
    def __init__(self):
        pass

    def _transform(self, sequence, class1, class2, class3):
        char = ''
        for element in sequence:
            if element in class1:
                char += 'P'
            elif element in class2:
                char += 'N'
            elif element in class3:
                char += 'H'
            else:
                char += 'X'
        return char

    def _computing_CTD(self, PNH_sequence):
        CTD_feature_vector = []

        C_P = PNH_sequence.count('P') / len(PNH_sequence)
        CTD_feature_vector.append(C_P)
        C_N = PNH_sequence.count('N') / len(PNH_sequence)
        CTD_feature_vector.append(C_N)
        C_H = PNH_sequence.count('H') / len(PNH_sequence)
        CTD_feature_vector.append(C_H)

        T_PN_NP = (PNH_sequence.count('PN') + PNH_sequence.count('NP')) / (len(PNH_sequence) - 1)
        CTD_feature_vector.append(T_PN_NP)
        T_PH_HP = (PNH_sequence.count('PH') + PNH_sequence.count('HP')) / (len(PNH_sequence) - 1)
        CTD_feature_vector.append(T_PH_HP)
        T_NH_HN = (PNH_sequence.count('NH') + PNH_sequence.count('HN')) / (len(PNH_sequence) - 1)
        CTD_feature_vector.append(T_NH_HN)

        P_N_H_LIST = [list() for i in range(3)]
        for i, v in enumerate(PNH_sequence):
            if v == 'P':
                P_N_H_LIST[0].append(i)
            elif v == 'N':
                P_N_H_LIST[1].append(i)
            elif v == 'H':
                P_N_H_LIST[2].append(i)

        P_N_H_per = [0 for i in range(15)]
        for i in range(3):
            if P_N_H_LIST[i]:
                P_N_H_per[i*5] = (P_N_H_LIST[i][0] + 1) / len(PNH_sequence)
                P_N_H_per[i*5 + 1] = (P_N_H_LIST[i][int((len(P_N_H_LIST[i]) * 0.25)) - 1] + 1) / len(PNH_sequence)
                P_N_H_per[i*5 + 2] = (P_N_H_LIST[i][int((len(P_N_H_LIST[i]) * 0.5)) - 1] + 1) / len(PNH_sequence)
                P_N_H_per[i*5 + 3] = (P_N_H_LIST[i][int((len(P_N_H_LIST[i]) * 0.75)) - 1] + 1) / len(PNH_sequence)
                P_N_H_per[i*5 + 4] = (P_N_H_LIST[i][-1] + 1) / len(PNH_sequence)
            else:
                for j in range(i*5, i*5+5):
                    P_N_H_per[j] = 0
        for i in range(15):
            CTD_feature_vector.append(P_N_H_per[i])

        # print(CTD_feature_vector[0])
        return CTD_feature_vector

    def get_feature_name(self):
        PNH_list = ['P', 'N', 'H']
        CTD_list = ['C', 'T', 'D']
        T_list = ['PN_NP', 'PH_HP', 'NH_HN']
        D_list = ['0', '25', '50', '75', '100']
        first_row_name = []
        for i in range(1, 8):
            for j in CTD_list:
                if j == 'C':
                    for k1 in PNH_list:
                        char = ('CTD_%a_' + j + '_' + k1) % i
                        first_row_name.append(char)
                elif j == 'T':
                    for k2 in T_list:
                        char = ('CTD_%a_' + j + '_' + k2) % i
                        first_row_name.append(char)
                elif j == 'D':
                    for k3 in PNH_list:
                        for l in D_list:
                            char = ('CTD_%a_' + j + '_' + k3 + '_' + l) % i
                            first_row_name.append(char)
        return first_row_name

    def main(self, sequence):
        feature1 = self._computing_CTD(self._transform(sequence,'RKEDQN','GASTPHY','CLVIMFW'))
        feature2 = self._computing_CTD(self._transform(sequence,'GASTPDC','NVEQIL','MHKFRYW'))
        feature3 = self._computing_CTD(self._transform(sequence,'LIFWCMVY','PATGS','HQRKNED'))
        feature4 = self._computing_CTD(self._transform(sequence,'GASDT','CPNVEQIL','KMHFRYW'))
        feature5 = self._computing_CTD(self._transform(sequence,'KR','ANCQGHILMFPSTWYV','DE'))
        feature6 = self._computing_CTD(self._transform(sequence,'EALMQKRH','VIYCWFT','GNPSD'))
        feature7 = self._computing_CTD(self._transform(sequence,'ALFCGIVW','PKQEND','MRSTHY'))
        feature_total = feature1 + feature2 + feature3 + feature4 + feature5 + feature6 + feature7
        feature = get_the_float(feature_total)
        # print(len(feature))
        return feature


class CC:
    def __init__(self):
        pass


    def append_l1_l2_l3(self, list1, list2, list3, num1, num2, num3):
        list1.append(num1)
        list2.append(num2)
        list3.append(num3)

    def get_feature_name(self):
        feature_name_list = ['Type1Xmean', 'Type1Ymean', 'Type1Zmean',
                             'Type2Xmean', 'Type2Ymean', 'Type2Zmean',
                             'Type3Xmean', 'Type3Ymean', 'Type3Zmean',
                             'Type4Xmean', 'Type4Ymean', 'Type4Zmean',
                             'AllXmean', 'AllYmean', 'AllZmean',
                             'CumXmean', 'CumYmean', 'CumZmean']
        return feature_name_list

    def main(self, sequence):
        lst = []
        seq_length = len(sequence)
        AA_char = 'AVLIPFWMGSTCYNQKRHDE'
        temp_num_list = [8, 8, 8, 7, 7, 7, 3, 3, 3, 2, 2, 2]
        Fs = []

        for char in AA_char:
            Fs.append(sequence.count(char))

        values = [0] * 12
        values[0] = Fs[0] * (-36.6794) + Fs[1] * (-66.6159) + Fs[2] * (-94.4033) + Fs[3] * (-100.0599) + Fs[4] * (-84.6835) + Fs[5] * (
            -101.0966) + Fs[6] * (-12.5874) + Fs[7] * (-111.3198)
        values[1] = Fs[0] * (58.8302) + Fs[1] * (62.1625) + Fs[2] * (38.8595) + Fs[3] * (20.2361) + Fs[4] * (29.1446) + Fs[5] * (
            -79.3826) + Fs[6] * (-158.3870) + Fs[7] * (-32.9407)
        values[2] = Fs[0] * (-55.9682) + Fs[1] * (-73.5564) + Fs[2] * (-82.4134) + Fs[3] * (-82.4134) + Fs[4] * (-72.3001) + Fs[5] * (
            -103.7705) + Fs[6] * (-128.2684) + Fs[7] * (-93.7201)
        values[3] = Fs[8] * (-24.7561) + Fs[9] * (-43.6432) + Fs[10] * (-25.2154) + Fs[11] * (-52.4159) + Fs[12] * (
            -53.3563) + Fs[13] * (-44.2471) + Fs[14] * (-65.6344)
        values[4] = Fs[8] * (26.2092) + Fs[9] * (25.3163) + Fs[10] * (51.3147) + Fs[11] * (-25.2563) + Fs[12] * (-68.7011) + Fs[13] * (
            -45.4289) + Fs[14] * (-24.8605)
        values[5] = Fs[8] * (65.8804) + Fs[9] * (92.1974) + Fs[10] * (104.4787) + Fs[11] * (106.3209) + Fs[12] * (
            158.9550) + Fs[13] * (115.8828) + Fs[14] * (128.2518)
        values[6] = Fs[15] * (-3.6804) + Fs[16] * (-2.1627) + Fs[17] * (-0.6273)
        values[7] = Fs[15] * (-1.8878) + Fs[16] * (-4.4286) + Fs[17] * (4.3459)
        values[8] = Fs[15] * (-146.1415) + Fs[16] * (-174.1303) + Fs[17] * (-155.1379)
        values[9] = Fs[18] * (-15.3220) + Fs[19] * (-16.9336)
        values[10] = Fs[18] * (43.3370) + Fs[19] * (-47.8954)
        values[11] = Fs[18] * (-124.9110) + Fs[19] * (-138.0496)
        values_temp = [i for i in values]
        for i in range(12):
            values[i] = values[i]/temp_num_list[i]
            lst.append(values[i])
        AllXmean = (values_temp[0] + values_temp[3] + values_temp[6] + values_temp[9]) / seq_length
        lst.append(AllXmean)
        AllYmean = (values_temp[1] + values_temp[4] + values_temp[7] + values_temp[10]) / seq_length
        lst.append(AllYmean)
        AllZmean = (values_temp[2] + values_temp[5] + values_temp[8] + values_temp[11]) / seq_length
        lst.append(AllZmean)

        lst1 = []
        lst2 = []
        lst3 = []
        num_char = [
            ['-36.6794','58.8302','-55.9682'],
            ['-66.6159', '62.1625', '73.5564'],
            ['-94.4003', '38.8595', '-82.4134'],
            ['-100.5090', '20.2361', '-82.4134'],
            ['-84.6835', '29.1446', '-72.3001'],
            ['-101.0996', '-79.3826', '-103.7705'],
            ['-12.5874', '-158.3870', '-128.2684'],
            ['-111.3198', '-32.9407', '-93.7201'],
            ['-24.7561', '26.2092', '65.8804'],
            ['-43.6432', '25.313163', '92.1974'],
            ['-25.2154', '51.3147', '104.3787'],
            ['-51.4159', '-25.2563', '106.3209'],
            ['-53.3563', '-68.7011', '158.9550'],
            ['-44.2471', '-45.4280', '115.8828'],
            ['-65.6344', '-24.8605', '128.2518'],
            ['-3.6804', '-1.8878', '-146.1415'],
            ['-2.1627', '-4.4286', '-174.1303'],
            ['-0.6273', '4.3459', '-155.1379'],
            ['-15.3220', '43.3370', '-124.9110'],
            ['-16.9336', '-47.8954', '-138.0496']
        ]
        for i in range(0, len(sequence)):
            for j in range(len(AA_char)):
                if sequence[i] == 'X':
                    self.append_l1_l2_l3(lst1, lst2, lst3, '0', '0', '0')
                    break
                if sequence[i] == AA_char[j]:
                    self.append_l1_l2_l3(lst1, lst2, lst3, num_char[j][0], num_char[j][1], num_char[j][2])
        cumulative_sum_x = 0
        cumulative_sum_y = 0
        cumulative_sum_z = 0
        feature_total = []
        feature_total.extend(lst)
        for inx, val in enumerate(lst1):
            cumulative_sum_x += float(val) * (len(lst1) - inx)
        feature_total.append(cumulative_sum_x / len(lst1))
        for inx, val in enumerate(lst2):
            cumulative_sum_y += float(val) * (len(lst2) - inx)
        feature_total.append(cumulative_sum_y / len(lst2))

        for inx, val in enumerate(lst3):
            cumulative_sum_z += float(val) * (len(lst3) - inx)
        feature_total.append(cumulative_sum_z / len(lst3))

        feature = get_the_float(feature_total)

        return feature

def main():
    f_positive = open(sys.argv[1], 'r')
    f_negative = open(sys.argv[2], 'r')

    ###############修改部分############
    # 获取长度
    length = 0
    for ele_pos in f_positive:
        if ele_pos[0] != '>':
            pos_sequence = ele_pos.strip()
            length = len(pos_sequence)
            break
    ###############修改部分############

    f_positive = open(sys.argv[1], 'r')
    f_negative = open(sys.argv[2], 'r')

    f_total = open('all_feature.csv', 'w')
    name_all_feature = []
    name_all_feature.extend(PCP9().get_feature_name(length))
    name_all_feature.extend(PWAA().get_feature_name())
    name_all_feature.extend(ARPC().get_feature_name(length))
    name_all_feature.extend(K_space(5).get_feature_name())
    name_all_feature.extend(EBGW(5).get_feature_name())
    name_all_feature.extend(CTD().get_feature_name())
    name_all_feature.extend(CC().get_feature_name())
    f_total.write('label,')
    for ix, name in enumerate(name_all_feature):
        f_total.write(name)

        if ix != (len(name_all_feature) - 1):
            f_total.write(',')
        else:

            f_total.write('\n')

    for ele_pos in f_positive:
        if ele_pos[0] != '>':
            f_total.write('1')
            f_total.write(',')
            feature_total = []
            pos_sequence = ele_pos.strip()
            feature_total.extend(PCP9().main(pos_sequence))
            feature_total.extend(PWAA().main(pos_sequence))
            feature_total.extend(ARPC().main(pos_sequence))
            for i in range(6):
                feature_total.extend(K_space(i).main(pos_sequence))
            for i in range(1,6):
                feature_total.extend(EBGW(i).main(pos_sequence))
            feature_total.extend(CTD().main(pos_sequence))

            feature_total.extend(CC().main(pos_sequence))

            for x, pos_feature in enumerate(feature_total):
                f_total.write(str(pos_feature))
                if x != (len(feature_total) - 1):
                    f_total.write(',')
                else:

                    f_total.write('\n')

    for ele_neg in f_negative:
        if ele_neg[0] != '>':

            f_total.write('0')
            f_total.write(',')
            feature_total = []
            neg_sequence = ele_neg.strip()
            feature_total.extend(PCP9().main(neg_sequence))
            feature_total.extend(PWAA().main(neg_sequence))
            feature_total.extend(ARPC().main(neg_sequence))
            for i in range(6):
                feature_total.extend(K_space(i).main(neg_sequence))
            for i in range(1,6):
                feature_total.extend(EBGW(i).main(neg_sequence))

            feature_total.extend(CTD().main(neg_sequence))
            feature_total.extend(CC().main(neg_sequence))
            for x, neg_feature in enumerate(feature_total):
                f_total.write(str(neg_feature))

                if x == (len(feature_total) - 1):
                    f_total.write('\n')
                else:
                    f_total.write(',')

if __name__ == '__main__':
    main()
