#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <random>
#include <sstream>
#include <iomanip>
#include "mini-gmp.h"
#include "mini-mpq.h"

#define KEY_SIZE 2048

using namespace std;

int gcd(int a, int b)
{
    int c;
    c = a % b;       //a>b
    while (c != 0)
    {
        a = b;
        b = c;
        c = a % b;
    }
    return b;
}


long long ext_gcd(long long a, long long b, long long& x, long long& y) {
    x = 1, y = 0;
    long long x2 = 0, y2 = 1;
    while (b) {
        long long q = a / b;
        long long r = a - q * b;
        a = b, b = r;

        long long tmp = x2;
        x2 = x - q * x2;
        x = tmp;

        tmp = y2;
        y2 = y - q * y2;
        y = tmp;
    }
 
    return a;
}


// 求 a 在模 m 意义下的乘法逆元
// 如果不存在乘法逆元，则返回 -1
long long mod_inv(long long a, long long m) {
    long long x, y;
    long long gcd_res = ext_gcd(a, m, x, y);
    if (gcd_res != 1) {
        // a 和 m 不互质，不存在乘法逆元
        return -1;
    }
    else {
        // 根据扩展欧几里得算法的结果计算乘法逆元
        return (x % m + m) % m;
    }
}

void mpz_mod_inv(mpz_t r, const mpz_t a, const mpz_t m)
{
    mpz_t x, y, gcd_res;
    mpz_init(x);
    mpz_init(y);
    mpz_init(gcd_res);

    mpz_gcdext(gcd_res, x, y,a,m);
    if (mpz_cmp_ui(gcd_res, 1) != 0) {
        // a 和 m 不互质，不存在乘法逆元
        mpz_init_set_si(r, -1);
    }
    else {
        // 根据扩展欧几里得算法的结果计算乘法逆元
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mod(tmp, x, m);
        mpz_add(r, tmp, m);
        mpz_mod(tmp, r, m);
        mpz_swap(tmp, r);
        mpz_clear(tmp);
        
    }
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(gcd_res);
}



// 计算幂余，即 x^y mod m
long long pow_mod(long long x, long long y, long long m)
{
    long long res = 1 % m;
    
    while (y > 0) {
        if (y % 2) res = (res * x) % m;
        y /= 2;
        x = (x * x) % m;
    }
    return res;
}



// 判断一个数是否为素数
bool is_prime(long long n)
{
    if (n <= 1) return false;
    for (long long i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

// 生成素数
long long generate_prime(long long min_val, long long max_val)
{
    long long p = 0;
    while (true) {
        p = rand() % (max_val - min_val + 1) + min_val;
        if (is_prime(p)) return p;
    }
}

void mpz_rand(mpz_t r, const mpz_t min_val, const mpz_t max_val)
{
    mpz_t interval,tmp_r;
    mpq_t rand_interval, rand_num, tmp, tmp2;
    mpz_init(interval);
    mpz_init(tmp_r);
    mpq_init(rand_interval);
    mpq_init(rand_num);
    mpq_init(tmp);
    mpq_init(tmp2);
    mpz_sub(interval, max_val, min_val);
    
    unsigned long rand_nums = 100; //增强随机性
    std::random_device e;
    std::uniform_real_distribution<double> u(0, 1); //随机数分布对象 
    std::minstd_rand linearRan(e());
    for (int i = 0; i < rand_nums; ++i)
    {
        mpq_set_d(rand_num, u(linearRan));
        mpq_set_num(rand_interval, interval);
        mpq_mul(tmp, rand_interval, rand_num);
        mpq_add(tmp, tmp, tmp);
    }
    mpq_set_ui(tmp2, rand_nums,1);
    mpq_div(tmp, tmp, tmp2);
 



    
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);
    mpq_get_den(den, tmp);
    mpq_get_num(num, tmp); 
    mpz_cdiv_q(r, num, den);//取整
    mpz_clear(num);
    mpz_clear(den);

   
    mpz_add(tmp_r, r, min_val);
    mpz_swap(tmp_r, r);

    mpz_clear(interval);
    mpz_clear(tmp_r);
    mpq_clear(rand_interval);
    mpq_clear(rand_num);
    mpq_clear(tmp);
    mpq_clear(tmp2);
    
     
    
    
}

void mpz_generate_prime(mpz_t p, const mpz_t min_val, const mpz_t max_val)
{
    mpz_init_set_ui(p, 0);
    mpz_t tmp;
    mpz_init(tmp);
    while (true) {
        mpz_rand(p, min_val, max_val);
        unsigned long s = mpz_mod_ui(tmp,p, 2);
        mpz_add_ui(p, p, s + 1); //保证是奇数
        int is_prime = mpz_probab_prime_p(p, 50);
        if (is_prime >0 ) break;
    }
    mpz_clear(tmp);
}

// 计算欧拉函数值
long long euler_func(long long p, long long q)
{
    return (p - 1) * (q - 1);
}

void mpz_euler_func(mpz_t r, const mpz_t p, const mpz_t q)
{
    mpz_t i, j, temp;
    mpz_init(i);                        //初始化i,j
    mpz_init(j);
    mpz_init_set_ui(temp,1);    //temp初始化赋值为1
    mpz_sub(i, p, temp);                //i = p - temp，即i=p-1
    mpz_sub(j, q, temp);                //j = q - temp，即j=q-1
    mpz_mul(r, j, i);              //*euler=j * i            
    mpz_clear(i);
    mpz_clear(j);
    mpz_clear(temp);

}

// 生成公钥和私钥
void generate_key(long long& e, long long& d, long long& n)
{
    // 随机生成两个素数 p 和 q
    srand(time(nullptr));
    long long p = generate_prime(1000, 100000);
    long long q = generate_prime(1000, 100000);

    // 计算 n 和 φ(n)
    n = p * q;
    long long phi_n = euler_func(p, q);

    // 选择一个整数 e，满足 1 < e < φ(n)，且 e 与 φ(n) 互质
    do {
        e = rand() % (phi_n - 2) + 2;
    } while (gcd(e, phi_n) != 1);

    // 计算 e 的模反元素 d
    d = mod_inv(e, phi_n);
}

// 生成公钥和私钥
void mpz_generate_key(mpz_t e, mpz_t d, mpz_t n, long key_size)
{
    mpz_t p, q, phi_n;
    mpz_t min_val, max_val;
    mpz_t tmp;
    mpz_init(p);
    mpz_init(q);
    mpz_init(phi_n);
    mpz_init(min_val);
    mpz_init(max_val);
    mpz_init_set_ui(tmp, 2);
    mpz_pow_ui(min_val, tmp, key_size / 2 -1);
    mpz_pow_ui(max_val, tmp, key_size / 2);
    mpz_generate_prime(p, min_val, max_val);
    mpz_generate_prime(q, min_val, max_val);

    mpz_euler_func(phi_n, p, q);
    mpz_mul(n, p, q);
    mpz_init_set_ui(e, 65537);

    // 计算 e 的模反元素 d
    mpz_mod_inv(d,e, phi_n);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(phi_n);
    mpz_clear(min_val);
    mpz_clear(max_val);
    mpz_clear(tmp);
}


// 加密
string encrypt(string msg, long long e, long long n)
{
    string res = "";
    for (char c : msg) {
        long long m = c;
        long long cip = pow_mod(m, e, n);
        res += to_string(cip) + " ";
    }
    return res;
}

//bytes to hex
string hex_encode(string str)
{
    std::stringstream ss;
 
    for (size_t i = 0; i < str.size(); ++i) {
        ss << std::hex << static_cast<int>(str[i]); // 将每个字符转换为16进制编码
    }

    std::string hex_str = ss.str(); // 获取16进制编码字符串
    return hex_str;
}

string hex_encode(vector<unsigned char> data, int offset, int length)
{
    std::stringstream ss;
    ss << std::hex << std::setfill('0');
    for (size_t i = 0; i < length; ++i) {
        ss << std::setw(2) << static_cast<int>(data[i + offset]); // 将每个字符转换为16进制编码
    }

    std::string hex_str = ss.str(); // 获取16进制编码字符串
    return hex_str;
}

//hex decode
string hex_decode(string hex_str)
{
    std::string str;
    if (hex_str.length() % 2 != 0)
    {
        hex_str.insert(hex_str.begin(), '0');
    }
    for (size_t i = 0; i < hex_str.length(); i += 2) {
        std::string byte_str = hex_str.substr(i, 2); // 获取两个字符作为一个字节的16进制编码字符串
        char byte = static_cast<char>(std::stoi(byte_str, nullptr, 16)); // 将16进制编码字符串转换为字节
        str += byte; // 将字节添加到字符串中
    }
    return str;
}

//hex decode
vector<unsigned char> hex_decode_bytes(string hex_str)
{
    vector<unsigned char> data;
    if (hex_str.length() % 2 != 0)
    {
        hex_str.insert(hex_str.begin(), '0');
    }
    for (size_t i = 0; i < hex_str.length(); i += 2) {
        std::string byte_str = hex_str.substr(i, 2); // 获取两个字符作为一个字节的16进制编码字符串
        unsigned char byte = static_cast<unsigned char>(std::stoi(byte_str, nullptr, 16)); // 将16进制编码字符串转换为字节
        data.push_back(byte);
    }
    return data;
}

// 加密
string mpz_encrypt(string msg, const mpz_t e, const mpz_t n)
{
    string res = "";
    mpz_t m, cip;
    mpz_init(m);
    mpz_init(cip);

    const int key_size = KEY_SIZE;
    const int block_size = key_size / 8 / 2;

    for (int i = 0; i < msg.length(); i+=block_size)
    {
        
        string hex_str;
        if (i + block_size > msg.length())
        {
            hex_str = hex_encode(msg.substr(i, msg.length() - i));
        }
        else
        {
           hex_str = hex_encode(msg.substr(i, block_size));
        }
        mpz_set_str(m, hex_str.c_str(), 16);
        mpz_powm(cip, m, e, n);
        char buff[KEY_SIZE * 2] = { 0 };
        mpz_get_str(buff, 16, cip);
        res += buff;
        res += " ";
       
    }

    mpz_clear(m);
    mpz_clear(cip);
    return res;
}


// 加密
string mpz_encrypt(vector<unsigned char> msg, const mpz_t e, const mpz_t n)
{
    string res = "";
    mpz_t m, cip;
    mpz_init(m);
    mpz_init(cip);

    const int key_size = KEY_SIZE;
    const int block_size = key_size / 8 / 2;

    for (int i = 0; i < msg.size(); i += block_size)
    {

        string hex_str;
        if (i + block_size > msg.size())
        {
            hex_str = hex_encode(msg, i, msg.size() - i);
        }
        else
        {
            hex_str = hex_encode(msg, i, block_size);
        }
        mpz_set_str(m, hex_str.c_str(), 16);
        mpz_powm(cip, m, e, n);
        char buff[KEY_SIZE * 2] = { 0 };
        mpz_get_str(buff, 16, cip);
        res += buff;
        res += " ";

    }

    mpz_clear(m);
    mpz_clear(cip);
    return res;
}


// 解密
string decrypt(string msg, long long d, long long n)
{
    string res = "";
    string num = "";
    for (char c : msg) {
        if (c == ' ') {
            long long cip = stoll(num);
            long long m = pow_mod(cip, d, n);
            res += (char)m;
            num = "";
        }
        else {
            num += c;
        }
    }
    return res;
}

// 解密
string mpz_decrypt(string msg, const mpz_t d, const mpz_t n)
{
    string res = "";
    string num = "";
    mpz_t m, cip;
    mpz_init(m);
    mpz_init(cip);
    for (char c : msg) {
        if (c == ' ') {
            mpz_init_set_str(cip, num.c_str(), 16);
            mpz_powm(m,cip, d, n);
            char buff[KEY_SIZE*2] = {0};
            mpz_get_str(buff,16,m);
            string origin_str = hex_decode(buff);
            res += origin_str;
            num = "";
        }
        else {
            num += c;
        }
    }
    return res;
}

// 解密
vector<unsigned char> mpz_decrypt_bytes(string msg, const mpz_t d, const mpz_t n)
{
    vector<unsigned char> res;
    string num = "";
    mpz_t m, cip;
    mpz_init(m);
    mpz_init(cip);
    for (char c : msg) {
        if (c == ' ') {
            mpz_init_set_str(cip, num.c_str(), 16);
            mpz_powm(m, cip, d, n);
            char buff[KEY_SIZE * 2] = { 0 };
            mpz_get_str(buff, 16, m);
            auto origin_str = hex_decode_bytes(buff);
            for (auto& t : origin_str)
            {
                res.push_back(t);
            }
            num = "";
        }
        else {
            num += c;
        }
    }
    return res;
}

int main()
{
    //long long e, d, n;
    //generate_key(e, d, n);
    //cout << "Public Key: " << e << ", " << n << endl;
    //cout << "Private Key: " << d << ", " << n << endl;

    //string msg = "Hello, world!";
    //cout << "Original Message: " << msg << endl;

    //string cipher_text = encrypt(msg, e, n);
    //cout << "Encrypted Message: " << cipher_text << endl;

    //string plain_text = decrypt(cipher_text, d, n);
    //cout << "Decrypted Message: " << plain_text << endl;



    mpz_t e, d, n;
    mpz_init(e);
    mpz_init(d);
    mpz_init(n);
    mpz_generate_key(e, d, n, KEY_SIZE);
    char buff[KEY_SIZE*2] = {0};
    cout << "n length" << strlen(mpz_get_str(buff, 16, n)) << endl;
    cout << "Public Key: " << mpz_get_str(buff, 16, e) << ", ";
    cout << mpz_get_str(buff, 16, n) <<endl;
    cout << "Private Key: length " <<strlen(mpz_get_str(buff, 16, d))<<" " << mpz_get_str(buff, 16, d) << ", ";
    cout << mpz_get_str(buff, 16, n) << endl;

    //string msg = "Hello, world!";
    //for (int i = 0; i < 2000; ++i)
    //{
    //    msg.push_back(rand() % 26 + 'A');
    //}
    //
    //cout << "Original Message: " << msg << endl;

    //string cipher_text = mpz_encrypt(msg, e, n);
    //cout << "Encrypted Message: " << cipher_text << endl;

    //string plain_text = mpz_decrypt(cipher_text, d, n);
    //cout << "Decrypted Message: " << plain_text << endl;

    //cout << (msg == plain_text) << endl;

    vector<unsigned char> data;
    for (int i=0;i<2000;++i)
    {
        data.push_back(rand()%256);
    }
    string encrypt_text = mpz_encrypt(data, e, n);
    auto decrypt_data = mpz_decrypt_bytes(encrypt_text, d, n);
    cout << "length check " << (data.size() == decrypt_data.size()) << endl;
    for (int i = 0; i < data.size(); ++i)
    {
        if (data[i] != decrypt_data[i])
        {
            cout << "fail " <<i << " " << (int)data[i] << " " << (int)decrypt_data[i] << endl;
            break;
        }
    }

    return 0;
}
   
