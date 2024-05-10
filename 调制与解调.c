//
//  main.c
//  选择页面
//
//  Created by 邹娜 on 2024/5/7.
//

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>

#define PI 3.14159265
#define SAMPLING_FREQUENCY 10000.0 // 采样频率
#define CARRIER_FREQUENCY 1000.0 // 载波频率
#define CARRIER_FREQUENCY_2 2000.0 // 第二个载波频率
#define MODULATION_INDEX 0.5 // 调制指数
#define SIGNAL_FREQUENCY 100.0 // 信号频率
#define SIGNAL_AMPLITUDE 1.0 // 信号幅度
#define BIT_RATE 100.0 // 比特率
#define MODULATION_ORDER 16 // QAM调制阶数
#define SYMBOLS_PER_BIT 2 // 每个比特对应的符号数


// 原始信号的频率
double originalFreq = 100.0;  // Hz

// 载波信号的频率
double carrierFreq = 1000.0;  // Hz

// 调制深度
double modulationIndex = 0.5;

// 采样率
double samplingRate = 10000.0;  // Hz

// 信号持续时间
double duration = 1.0;  // seconds
double t;


// 生成原始信号
void generate_signal(double *signal, int signal_length) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    for (int i = 0; i < signal_length; i++) {
        signal[i] = SIGNAL_AMPLITUDE * sin(2 * PI * SIGNAL_FREQUENCY * i * time_period);
    }
}

// 方波的占空比
double dutyCycle = 0.5;

// 方波信号函数
double squareWave(double t, double freq, double duty) {
    double period = 1.0 / freq;
    double t_mod = fmod(t, period);
    if (t_mod < period * duty) {
        return 1.0;
    } else {
        return -1.0;
    }
}

// 载波信号函数
double carrierSignal(double t) {
    return sin(2 * PI * carrierFreq * t);
}

// AM调制
void am_modulation(double *signal, int signal_length, double *modulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    for (int i = 0; i < signal_length; i++) {
        modulated_signal[i] = (1 + MODULATION_INDEX * signal[i]) * cos(2 * PI * CARRIER_FREQUENCY * i * time_period);
    }
}

// AM解调
void am_demodulation(double *modulated_signal, int signal_length, double *demodulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    for (int i = 0; i < signal_length; i++) {
        demodulated_signal[i] = modulated_signal[i] * cos(2 * PI * CARRIER_FREQUENCY * i * time_period);
    }
}

// FM调制
void fm_modulation(double *signal, int signal_length, double *modulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    double phase = 0.0; // 初始相位
    for (int i = 0; i < signal_length; i++) {
        phase += 2 * PI * CARRIER_FREQUENCY * time_period + MODULATION_INDEX * signal[i];
        modulated_signal[i] = cos(phase);
    }
}

// FM解调
void fm_demodulation(double *modulated_signal, int signal_length, double *demodulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    double previous_phase = 0.0;
    for (int i = 0; i < signal_length; i++) {
        double phase_difference = atan2(sin(modulated_signal[i] - previous_phase), cos(modulated_signal[i] - previous_phase));
        // 通过相位差计算解调信号
        demodulated_signal[i] = phase_difference / (2 * PI * time_period * MODULATION_INDEX);
        previous_phase = modulated_signal[i];
    }
}

// PM调制
void pm_modulation(double *signal, int signal_length, double *modulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    double phase = 0.0; // 初始相位
    for (int i = 0; i < signal_length; i++) {
        modulated_signal[i] = cos(2 * PI * CARRIER_FREQUENCY * i * time_period + MODULATION_INDEX * signal[i]);
    }
}

// PM解调
void pm_demodulation(double *modulated_signal, int signal_length, double *demodulated_signal) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    double previous_phase = 0.0;
    for (int i = 0; i < signal_length; i++) {
        double phase_difference = atan2(sin(2 * PI * CARRIER_FREQUENCY * i * time_period), cos(2 * PI * CARRIER_FREQUENCY * i * time_period)) - previous_phase;
        // 通过相位差计算解调信号
        demodulated_signal[i] = phase_difference / (2 * PI * time_period * MODULATION_INDEX);
        previous_phase = atan2(sin(2 * PI * CARRIER_FREQUENCY * i * time_period), cos(2 * PI * CARRIER_FREQUENCY * i * time_period));
    }
}
// 生成二进制数据流
void generate_binary_data(int *binary_data, int bit_count) {
    // 简单生成一个固定的二进制数据序列
    for (int i = 0; i < bit_count; i++) {
        binary_data[i] = rand() % 2; // 随机生成0或1
    }
}

// ASK调制
void ask_modulation(int *binary_data, int bit_count, double *modulated_signal, int signal_length) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        double amplitude = (binary_data[i] == 1) ? 1.0 : 0.0; // 根据比特值选择振幅
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                modulated_signal[index] = amplitude * sin(2 * PI * CARRIER_FREQUENCY * index * time_period);
            }
        }
    }
}

// ASK解调
void ask_demodulation(double *modulated_signal, int signal_length, int *demodulated_bits, int bit_count) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        double average_amplitude = 0.0;
        // 计算一个比特对应的信号的平均振幅
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                average_amplitude += fabs(modulated_signal[index]);
            }
        }
        average_amplitude /= samples_per_bit;

        // 如果平均振幅超过阈值，则认为该比特为1
        demodulated_bits[i] = (average_amplitude > 0.5) ? 1 : 0;
    }
}

// FSK调制
void fsk_modulation(int *binary_data, int bit_count, double *modulated_signal, int signal_length) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        double frequency = (binary_data[i] == 1) ? CARRIER_FREQUENCY : CARRIER_FREQUENCY_2; // 根据比特值选择载波频率
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                modulated_signal[index] = sin(2 * PI * frequency * index * time_period);
            }
        }
    }
}

// FSK解调
void fsk_demodulation(double *modulated_signal, int signal_length, int *demodulated_bits, int bit_count) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        // 计算每个比特的平均幅度
        double sum = 0.0;
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                sum += modulated_signal[index];
            }
        }
        double average_amplitude = sum / samples_per_bit;

        // 根据平均幅度判断解调后的比特值
        demodulated_bits[i] = (average_amplitude > 0) ? 1 : 0;
    }
}

// PSK调制
void psk_modulation(int *binary_data, int bit_count, double *modulated_signal, int signal_length) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        double phase = binary_data[i] ? PI : 0; // 根据比特值选择相位
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                modulated_signal[index] = sin(2 * PI * CARRIER_FREQUENCY * index * time_period + phase);
            }
        }
    }
}

// PSK解调
void psk_demodulation(double *modulated_signal, int signal_length, int *demodulated_bits, int bit_count) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    for (int i = 0; i < bit_count; i++) {
        // 计算每个比特的相位
        double sum_cos = 0.0;
        double sum_sin = 0.0;
        for (int j = 0; j < samples_per_bit; j++) {
            int index = i * samples_per_bit + j;
            if (index < signal_length) {
                sum_cos += modulated_signal[index] * cos(2 * PI * CARRIER_FREQUENCY * index * time_period);
                sum_sin += modulated_signal[index] * sin(2 * PI * CARRIER_FREQUENCY * index * time_period);
            }
        }
        double phase = atan2(sum_sin, sum_cos);

        // 根据相位判断解调后的比特值
        demodulated_bits[i] = (phase > PI / 2) ? 1 : 0;
    }
}

// QAM调制
void qam_modulation(int *binary_data, int bit_count, double *modulated_signal, int signal_length) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_symbol = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    int symbols_count = bit_count / SYMBOLS_PER_BIT;
    double constellation[MODULATION_ORDER][2]; // 符号星座
    // 初始化16-QAM符号星座
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            constellation[i * 4 + j][0] = (2 * i - 3) / sqrt(10);
            constellation[i * 4 + j][1] = (2 * j - 3) / sqrt(10);
        }
    }

    for (int i = 0; i < symbols_count; i++) {
        int index = 0;
        for (int j = 0; j < SYMBOLS_PER_BIT; j++) {
            index += binary_data[i * SYMBOLS_PER_BIT + j] << (SYMBOLS_PER_BIT - 1 - j);
        }
        double I = constellation[index][0];
        double Q = constellation[index][1];

        for (int j = 0; j < samples_per_symbol; j++) {
            int sample_index = i * samples_per_symbol + j;
            if (sample_index < signal_length) {
                double phase = 2 * PI * CARRIER_FREQUENCY * sample_index * time_period;
                modulated_signal[sample_index] = I * cos(phase) - Q * sin(phase);
            }
        }
    }
}

// QAM解调
void qam_demodulation(double *modulated_signal, int signal_length, int *demodulated_bits, int bit_count) {
    double time_period = 1.0 / SAMPLING_FREQUENCY;
    int samples_per_symbol = (int)(SAMPLING_FREQUENCY / BIT_RATE);

    int symbols_count = bit_count / SYMBOLS_PER_BIT;
    double constellation[MODULATION_ORDER][2]; // 符号星座
    // 初始化16-QAM符号星座
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            constellation[i * 4 + j][0] = (2 * i - 3) / sqrt(10);
            constellation[i * 4 + j][1] = (2 * j - 3) / sqrt(10);
        }
    }

    for (int i = 0; i < symbols_count; i++) {
        double max_distance = INFINITY;
        int detected_symbol = 0;

        for (int k = 0; k < MODULATION_ORDER; k++) {
            double I = constellation[k][0];
            double Q = constellation[k][1];

            double sum = 0.0;
            for (int j = 0; j < samples_per_symbol; j++) {
                int sample_index = i * samples_per_symbol + j;
                if (sample_index < signal_length) {
                    double phase = 2 * PI * CARRIER_FREQUENCY * sample_index * time_period;
                    sum += modulated_signal[sample_index] * (I * cos(phase) - Q * sin(phase));
                }
            }

            if (sum < max_distance) {
                max_distance = sum;
                detected_symbol = k;
            }
        }

        // 根据检测到的符号确定解调后的比特值
        for (int j = 0; j < SYMBOLS_PER_BIT; j++) {
            demodulated_bits[i * SYMBOLS_PER_BIT + j] = (detected_symbol >> (SYMBOLS_PER_BIT - 1 - j)) & 1;
        }
    }
}

// 保存信号数据到文件
void save_signal_to_file(double *signal, int signal_length, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file.\n");
        exit(1);
    }

    for (int i = 0; i < signal_length; i++) {
        fprintf(file, "%f\n", signal[i]);
    }
    
    fclose(file);
    
}
// 保存二进制数据到文件
void save_binary_data_to_file(int *binary_data, int bit_count, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file.\n");
        exit(1);
    }

    for (int i = 0; i < bit_count; i++) {
        fprintf(file, "%d\n", binary_data[i]);
    }

    fclose(file);
}

int main(int argc, const char * argv[]) {
    
    int choice;
    void AMoutput();
    void FMoutput();
    void PMoutput();
    void ASKoutput();
    void FSKoutput();
    void PSKoutput();
    void QAMoutput();
    do {
        // 显示选择页面
        printf("\n请选择调制方式：\n");
        printf("1. AM\n");
        printf("2. FM\n");
        printf("3. PM\n");
        printf("4. ASK\n");
        printf("5. FSK\n");
        printf("6. PSK\n");
        printf("7. QAM\n");
        printf("请输入选项号码：");

        // 获取用户输入的选项
        scanf("%d", &choice);

        // 根据用户的选择执行相应操作
        switch (choice) {
                    case 1:
                AMoutput();
                        break;
                    case 2:
                FMoutput();
                        break;
                    case 3:
                PMoutput();
                        break;
                    case 4:
                ASKoutput();
                        break;
                    case 5:
                FSKoutput();
                        break;
                    case 6:
                PSKoutput();
                        break;
                    case 7:
                QAMoutput();
                        break;
                    default:
                        printf("无效的选项，请重新选择。\n");
                        break;
                }


    } while (choice != 4);

    return 0;
}

void AMoutput(){
    int signal_length = 1000; // 信号长度
        double signal[signal_length]; // 原始信号
        double modulated_signal[signal_length]; // 调制信号
        double demodulated_signal[signal_length]; // 解调信号

    // 生成原始信号
        generate_signal(signal, signal_length);

        // 保存原始信号到文件
        save_signal_to_file(signal, signal_length, "original_signal.txt");

        // 调制
        am_modulation(signal, signal_length, modulated_signal);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 解调
        am_demodulation(modulated_signal, signal_length, demodulated_signal);

        // 保存解调信号到文件
        save_signal_to_file(demodulated_signal, signal_length, "demodulated_signal.txt");

        printf("AM modulation and demodulation completed.\n");

}

void FMoutput(){
    int signal_length = 1000; // 信号长度
        double signal[signal_length]; // 原始信号
        double modulated_signal[signal_length]; // 调制信号
        double demodulated_signal[signal_length]; // 解调信号

        // 生成原始信号
        generate_signal(signal, signal_length);

        // 保存原始信号到文件
        save_signal_to_file(signal, signal_length, "original_signal.txt");

        // 调制
        fm_modulation(signal, signal_length, modulated_signal);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 解调
        fm_demodulation(modulated_signal, signal_length, demodulated_signal);

        // 保存解调信号到文件
        save_signal_to_file(demodulated_signal, signal_length, "demodulated_signal.txt");

        printf("FM modulation and demodulation completed.\n");

}

void PMoutput(){
    int signal_length = 1000; // 信号长度
        double signal[signal_length]; // 原始信号
        double modulated_signal[signal_length]; // 调制信号
        double demodulated_signal[signal_length]; // 解调信号

        // 生成原始信号
        generate_signal(signal, signal_length);

        // 保存原始信号到文件
        save_signal_to_file(signal, signal_length, "original_signal.txt");

        // 调制
        pm_modulation(signal, signal_length, modulated_signal);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 解调
        pm_demodulation(modulated_signal, signal_length, demodulated_signal);

        // 保存解调信号到文件
        save_signal_to_file(demodulated_signal, signal_length, "demodulated_signal.txt");

        printf("PM modulation and demodulation completed.\n");

}

void ASKoutput(){
    int bit_count = 100; // 二进制数据的位数
        int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);
        int signal_length = bit_count * samples_per_bit; // 计算调制信号的总长度

        int binary_data[bit_count]; // 二进制数据
        double modulated_signal[signal_length]; // 调制信号
        int demodulated_bits[bit_count]; // 解调后的二进制数据

        // 生成随机二进制数据
        generate_binary_data(binary_data, bit_count);

        // 保存原始二进制数据到文件
        save_binary_data_to_file(binary_data, bit_count, "original_binary_data.txt");

        // 进行ASK调制
        ask_modulation(binary_data, bit_count, modulated_signal, signal_length);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 进行ASK解调
        ask_demodulation(modulated_signal, signal_length, demodulated_bits, bit_count);

        // 保存解调后的二进制数据到文件
        save_binary_data_to_file(demodulated_bits, bit_count, "demodulated_binary_data.txt");

        // 检查解调结果与原始数据是否一致
        int error_count = 0;
        for (int i = 0; i < bit_count; i++) {
            if (binary_data[i] != demodulated_bits[i]) {
                error_count++;
             }
        }
    printf("ASK modulation and demodulation completed.\n");
    printf("Bit errors: %d\n", error_count);
}

void FSKoutput(){
    int bit_count = 100; // 二进制数据的位数
        int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);
        int signal_length = bit_count * samples_per_bit; // 计算调制信号的总长度

        int binary_data[bit_count]; // 二进制数据
        double modulated_signal[signal_length]; // 调制信号
        int demodulated_bits[bit_count]; // 解调后的二进制数据

        // 生成随机二进制数据
        generate_binary_data(binary_data, bit_count);

        // 保存原始二进制数据到文件
        save_binary_data_to_file(binary_data, bit_count, "original_binary_data.txt");

        // 进行FSK调制
        fsk_modulation(binary_data, bit_count, modulated_signal, signal_length);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 进行FSK解调
        fsk_demodulation(modulated_signal, signal_length, demodulated_bits, bit_count);

        // 保存解调后的二进制数据到文件
        save_binary_data_to_file(demodulated_bits, bit_count, "demodulated_binary_data.txt");
    
        // 检查解调结果与原始数据是否一致
        int error_count = 0;
        for (int i = 0; i < bit_count; i++) {
            if (binary_data[i] != demodulated_bits[i]) {
                 error_count++;
             }
         }
    printf("FSK modulation and demodulation completed.\n");
    printf("Bit errors: %d\n", error_count);
}

void PSKoutput(){
    int bit_count = 100; // 二进制数据的位数
        int samples_per_bit = (int)(SAMPLING_FREQUENCY / BIT_RATE);
        int signal_length = bit_count * samples_per_bit; // 计算调制信号的总长度

        int binary_data[bit_count]; // 二进制数据
        double modulated_signal[signal_length]; // 调制信号
        int demodulated_bits[bit_count]; // 解调后的二进制数据

        // 生成随机二进制数据
        generate_binary_data(binary_data, bit_count);

        // 保存原始二进制数据到文件
        save_binary_data_to_file(binary_data, bit_count, "original_binary_data.txt");

        // 进行PSK调制
        psk_modulation(binary_data, bit_count, modulated_signal, signal_length);

        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 进行PSK解调
        psk_demodulation(modulated_signal, signal_length, demodulated_bits, bit_count);

        // 保存解调后的二进制数据到文件
        save_binary_data_to_file(demodulated_bits, bit_count, "demodulated_binary_data.txt");
       // 检查解调结果与原始数据是否一致
       int error_count = 0;
       for (int i = 0; i < bit_count; i++) {
           if (binary_data[i] != demodulated_bits[i]) {
               error_count++;
             }
         }
    
    printf("PSK modulation and demodulation completed.\n");
    printf("Bit errors: %d\n", error_count);

}

void QAMoutput(){
    int bit_count = 100; // 二进制数据的位数
        int samples_per_symbol = (int)(SAMPLING_FREQUENCY / BIT_RATE);
        int signal_length = bit_count * samples_per_symbol / SYMBOLS_PER_BIT; // 计算调制信号的总长度

        int binary_data[bit_count]; // 二进制数据
        double modulated_signal[signal_length]; // 调制信号
        int demodulated_bits[bit_count]; // 解调后的二进制数据

        // 生成随机二进制数据
        generate_binary_data(binary_data, bit_count);
    
        // 保存原始二进制数据到文件
        save_binary_data_to_file(binary_data, bit_count, "original_binary_data.txt");

        // 进行QAM调制
        qam_modulation(binary_data, bit_count, modulated_signal, signal_length);
    
        // 保存调制信号到文件
        save_signal_to_file(modulated_signal, signal_length, "modulated_signal.txt");

        // 进行QAM解调
        qam_demodulation(modulated_signal, signal_length, demodulated_bits, bit_count);

        // 保存解调后的二进制数据到文件
        save_binary_data_to_file(demodulated_bits, bit_count, "demodulated_binary_data.txt");
        // 检查解调结果与原始数据是否一致
        int error_count = 0;
        for (int i = 0; i < bit_count; i++) {
            if (binary_data[i] != demodulated_bits[i]) {
                error_count++;
            }
        }

    printf("QAM modulation and demodulation completed.\n");
    printf("Bit errors: %d\n", error_count);

}
