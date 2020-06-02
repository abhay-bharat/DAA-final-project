#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "intal.h"
#define SIZE 1000

// Library "intal" - Integer of arbitray length
// intal is a nonnegative integer of arbitrary length (not more than 1000 decimal digits).
// The integer is stored as a null terminated string of ASCII characters.
// String of decimal digits ('0' thru '9') are stored in big endian style.
// That is, the most significant digit is at the head of the string.
// Eg: Integer 25 is stored in str as '2' at str[0], '5' at str[1], and null char at str[2].

// To generate your own testcases, you can use https://www.omnicalculator.com/math
static int find_length(const char* s)
{
    int len = 0;
    while(s[len] != '\0')
        ++len;
    return len;
}

//Helper function to reverse the string
static void reverse_string(char *dest,const char *src)
{
    int len = strlen(src);
    int i = 0;
    while(i != len)
    {
        dest[i] = src[len-i-1];
        ++i;
    }
    dest[len] = '\0';
}
 //Swapping the pointers
static void swap(char **a, char **b)
{
    char *temp = *a;
    *a = *b;
    *b = temp;
}
// Returns the sum of two intals.
char* intal_add(const char* intal1, const char* intal2)
{
    int n1, n2;
    int res_len = 0; //Length of resulting intal
    int d1, d2;
    char *temp_result;
    char *result;
    n1 = find_length(intal1);
    n2 = find_length(intal2);
    res_len = (n1 > n2)?n1:n2;
    int len = res_len;
    temp_result = (char*)malloc(res_len*sizeof(char)); 
    //Implement addition here
    int carry = 0;
    int ans = 0;
    while(n1 != 0 && n2 != 0)
    {
        d1 = intal1[n1-1] - '0';
        d2 = intal2[n2-1] - '0';
        ans = d1 + d2 + carry;
        carry = ans/10;
        temp_result[res_len-1] = ans%10 + '0';
        --res_len;
        --n1;
        --n2;
    }
    if(n1 == 0 && n2 != 0)
    {
        while(n2 != 0)
        {           
            d2 = intal2[n2-1] - '0';
            ans = d2 + carry;
            carry = ans/10;
            temp_result[res_len-1] = ans%10 + '0';
            --res_len;
            --n2;
        }
    }
    else if(n2 == 0 && n1 != 0)
    {
        while(n1 != 0)
        {
            d1 = intal1[n1-1] - '0';
            ans = d1 + carry;
            carry = ans/10;
            temp_result[res_len-1] = ans%10 + '0';
            --res_len;
            --n1;
        }
    }
    if(carry > 0)
    {
        result = (char*)malloc((len+2)*sizeof(char));
        result[0] = carry + '0';
        int i = 1;
        while(i <= len)
        {
            result[i] = temp_result[i-1];
            ++i;
        }
        result[len+1] = '\0';
    }
    else
    {
        result = (char*)malloc((len+1)*sizeof(char));
        int i = 0;
        while(i < len)
        {
            result[i] = temp_result[i];
            ++i;
        }
        result[len] = '\0';
    }
    free(temp_result);
    return result;
}

// Returns the comparison value of two intals.
// Returns 0 when both are equal.
// Returns +1 when intal1 is greater, and -1 when intal2 is greater.
int intal_compare(const char* intal1, const char* intal2)
{
    int n1 = find_length(intal1);
    int n2 = find_length(intal2);
    if(n1 > n2)
        return 1;
    else if(n1 < n2)
        return -1;
    else if(n1 == n2)
    {
        for(int i =0;i<n1;++i)
        {
            if(intal1[i] > intal2[i])
                return 1;
            else if(intal1[i] < intal2[i])
                return -1;
        }
    }
     return 0; //When completes for loop and both conditions inside for loop failed, which means they are equal
}


// Returns the difference (obviously, nonnegative) of two intals.
char* intal_diff(const char* intal1, const char* intal2)
{
    int comp = intal_compare(intal1, intal2);
    char *res;
    char *s1, *s2;
    int n1,n2;
    n1 = find_length(intal1);
    n2 = find_length(intal2);
    if(comp == 0)
    {
        res = (char*)malloc(2*sizeof(char));
        strcpy(res,"0");
        return res;
    }
    else if(comp == -1)
    {
        s1 = (char*)malloc((n2+1)*sizeof(char));
        s2 = (char*)malloc((n1+1)*sizeof(char));
        reverse_string(s1, intal2);  //reverse string intal2 and store it in s1
        reverse_string(s2, intal1);
    }
    else if(comp == 1)
    {
        s1 = (char*)malloc((n1+1)*sizeof(char));
        s2 = (char*)malloc((n2+1)*sizeof(char));
        reverse_string(s1, intal1);
        reverse_string(s2, intal2);
    }
    int len = (n1>n2)?n1:n2;
    res = (char*)malloc((len+1)*sizeof(char));
    int count = 0;
    int sub;
    int carry = 0;
    for(int i = 0; i<(n1+n2-len);++i)
    {
        sub = (s1[i]-'0') - (s2[i]-'0') - carry;
        if(sub < 0)
        {
            sub = sub + 10;
            carry = 1;
        }
        else
            carry = 0;
        res[count] = sub + '0';
        count++;
    }
    for(int i = (n1+n2-len); i< len;++i)
    {
        sub = (s1[i]-'0') - carry;
        if(sub < 0)
        {
            sub = sub+10;
            carry = 1;
        }
        else
            carry = 0;
        res[count] = sub + '0';
        count++;
    }
    res[count] = '\0';
    reverse_string(s1, res); //Now s1 will hold the result
    free(s2);
    //Now store the final result without leading zeroes in the array res
    int i = 0;
    while(s1[i] == '0')
        ++i;
    int j = 0;
    while(s1[i] != '\0')
    {
        res[j] = s1[i];
        ++i;
        ++j;
    }
    res[j] = '\0';
    free(s1);
    return res;
}

// Returns the product of two intals.
char* intal_multiply(const char* intal1, const char* intal2)
{
    int n1 = find_length(intal1);
    int n2 = find_length(intal2);
    int *res = (int*)calloc((n1+n2),sizeof(int)); //This is to store answer in reverse order
    char *result = (char*)malloc((n1+n2+1)*sizeof(char));
    int l1 = 0;
    int l2 = 0;

    for (int i=n1-1; i>=0; i--)
    {
        int carry = 0;
        int d1 = intal1[i] - '0';

        l2 = 0;
        int sum = 0;
        for (int j= n2-1; j>=0; j--)
        {
            // Take current digit of second number
            int d2 = intal2[j] - '0';

            //Multipling current digit of first number with all digits of 2nd number and storing result by adding with previous result
            sum = d1*d2 + res[l1 + l2] + carry;

            // Carry for next iteration
            carry = sum/10;

            // Store result
            res[l1 + l2] = sum % 10;

            l2++;
        }
        if (carry > 0)
            res[l1 + l2] += carry;

        l1++; //shifting position of intal1 by 1 digit
   }

    // ignore '0's from the right
    int i = n1+n2 -1;
    while (i>=0 && res[i] == 0)
        i--;

    //When all numbers are 0 then return string "0"
    if (i == -1)
    {
        strcpy(result,"0");
        return result;
    }
    //Now store the final answer from the reverse ordered integer array into the string array in proper order
    int j = 0;
    while(i >= 0)
    {
        result[j] = res[i] + '0';
        --i;
        ++j;
    }
    result[j] = '\0';
    free(res);
    return result;
}
//helper function to compute the intermediate modulus using subtraction method
static void compute_inter_mod(char *intal1, const char *intal2, char *ans)
{
    int cmp = intal_compare(intal1, intal2);
    char *temp = NULL;
    while(cmp > 0)
    {
        temp = intal_diff(intal1, intal2);
        strcpy(intal1, temp);
        free(temp);
        cmp = intal_compare(intal1, intal2);
    }
    if(cmp == 0)
    {
        strcpy(ans,"0");
    }
    else
    {
        strcpy(ans,intal1);
    }
}

// Returns intal1 mod intal2
// The mod value should be in the range [0, intal2 - 1].
// intal2 > 1
// Implement a O(log intal1) time taking algorithm.
// O(intal1 / intal2) time taking algorithm may exceed time limit.
// O(intal1 / intal2) algorithm may repeatedly subtract intal2 from intal1.
// That will take intal1/intal2 iterations.
// You need to design a O(log intal1) time taking algorithm.
// Generate your own testcases at https://www.omnicalculator.com/math/modulo
char* intal_mod(const char* intal1, const char* intal2)
{
    char *res = (char*)malloc((SIZE+1)*sizeof(char));
    char *temp = (char*)malloc((SIZE+1)*sizeof(char));
    char *temp_intal1 = (char*)malloc((SIZE+1)*sizeof(char));
    int compare;
    compare = intal_compare(intal1, intal2);
    if(compare == -1)
    {
        strcpy(res, intal1);
        free(temp);
        free(temp_intal1);
        return res;
    }
    if(compare == 0)
    {
        strcpy(res, "0");
        free(temp);
        free(temp_intal1);
        return res;
    }
    else
    {
        int n = find_length(intal2);
        strcpy(temp_intal1, intal1);
        int n2 = n;
        while(compare > 0)
        {
            strncpy(temp, temp_intal1, n2); //Copy only certain size from starting of the string
            temp[n2] = '\0';
            if(intal_compare(intal2, temp) >= 0)
            {
                n2++;
            }
            else
            {
                compute_inter_mod(temp, intal2,res);

                strcpy(temp, temp_intal1+n2);
                strcpy(temp_intal1, temp);
                strcpy(temp, res); //Let temp hold answer 
                if(strcmp(res, "0") != 0)
                {
                    strcat(temp, temp_intal1); //concatinating res to beginning of temp_intal1
                    strcpy(temp_intal1, temp); //Copy back it to temp_intal1
                }                   
                compare = intal_compare(temp_intal1, intal2);
                n2 = n;
            }
        }
        if(compare == 0 || temp_intal1[0] == '\0')
        {
            strcpy(temp_intal1, "0");
        }
    }
    free(res);
    free(temp);
    return temp_intal1;
}

// Returns intal1 ^ intal2.
// Let 0 ^ n = 0, where n is an intal.
// Implement a O(log n) intal multiplications algorithm.
// 2^3000 has less than 1000 decimal digits. 3000 intal multiplications may exceed time limit.
char* intal_pow(const char* intal1, unsigned int n)
{
    char *res = (char*)malloc((SIZE*10)*sizeof(char));
    if(strcmp(intal1,"0")== 0)
    {
        strcpy(res, "0");
        return res;
    }
    if(n == 0)
    {
        strcpy(res, "1");
        return res;
    }
    char *temp = NULL;
    strcpy(res ,"1");
    for(int i =0;i< n/2;++i)
    {
        temp = intal_multiply(res, intal1);
        strcpy(res,temp);
        free(temp);
    }
    temp = intal_multiply(res, res);
    strcpy(res,temp);
    free(temp);
    if(n%2 != 0)
    {
        temp = intal_multiply(res, intal1);
        strcpy(res,temp);
        free(temp);
    }
    return res;
}

static char* comp_gcd(char *a, char *b)
{
    if(strcmp(b, "0") == 0)
    {
        char *res = (char*)malloc((strlen(a)+1)*sizeof(char));
        strcpy(res, a);
        free(a);
        free(b);
        return res;
    }
    char *temp = intal_mod(a,b);
    free(a);
    return comp_gcd(b, temp);
}
// Returns Greatest Common Devisor of intal1 and intal2.
// Let GCD be "0" if both intal1 and intal2 are "0" even though it is undefined, mathematically.
// Use Euclid's theorem to not exceed the time limit.
char* intal_gcd(const char* intal1, const char* intal2)
{
    if(strcmp(intal1, "0")==0 && strcmp(intal2, "0")==0)
    {
        char *res = (char*)malloc(2*sizeof(char));
        strcpy(res, "0");
        return res;
    }
    if(strcmp(intal1, "1") == 0 || strcmp(intal2, "1") == 0)
    {
        char *res = (char*)malloc(2*sizeof(char));
        strcpy(res, "1");
        return res;
    }

    char *temp_intal1 = (char*)malloc((strlen(intal1)+1)*sizeof(char));
    char *temp_intal2 = (char*)malloc((strlen(intal2)+1)*sizeof(char));

    strcpy(temp_intal1, intal1);
    strcpy(temp_intal2, intal2);

    return comp_gcd(temp_intal1, temp_intal2);
}

// Returns nth fibonacci number.
// intal_fibonacci(0) = intal "0".
// intal_fibonacci(1) = intal "1".
char* intal_fibonacci(unsigned int n)
{
    char *res = (char*)malloc((SIZE*10)*sizeof(char));
    if(n == 0)
    {
        strcpy(res, "0");
        return res;
    }
    if(n == 1)
    {
        strcpy(res, "1");
        return res;
    }
    char *temp = NULL;

    char *prev1 = (char*)malloc((SIZE*10)*sizeof(char));
    strcpy(prev1, "0");
    char *prev2 = (char*)malloc((SIZE*10)*sizeof(char));
    strcpy(prev2, "1");
    for(int i = 0;i <= n-2;++i)
    {
        temp = intal_add(prev1, prev2);
        strcpy(res, temp);
        free(temp);
        strcpy(prev1, prev2);
        strcpy(prev2, res);
    }
    free(prev1);
    free(prev2);
    return res;
}

static void num_to_string(char *str, unsigned int n)
{
    int i =0;
    while(n != 0)
    {
        str[i] = (n%10) + '0';
        n = n/10;
        ++i; 
    }
    str[i] = '\0';
}

// Returns the factorial of n.
char* intal_factorial(unsigned int n)
{
    char *res = (char*)malloc((SIZE*10)*sizeof(char));
    char *temp = NULL;
    if(n == 0)
    {
        strcpy(res,"1");
        return res;
    }

    char rev_num[SIZE*10];  //To hold the reverse of a number
    char num[SIZE*10];

    strcpy(res,"1");
    for(int i = 1; i<= n;++i)
    {
        //Converting number to string of digits
        num_to_string(rev_num, i);
        //Reversing the string to get original number in intal format
        reverse_string(num, rev_num); 
        temp = intal_multiply(num ,res); 
        strcpy(res,temp);
        free(temp);
    }
    return res;   
}

// Returns the Binomial Coefficient C(n,k).
// 0 <= k <= n
// C(n,k) < 10^1000 because the returning value is expected to be less than 10^1000.
// Use the Pascal's identity C(n,k) = C(n-1,k) + C(n-1,k-1)
// Make sure the intermediate intal values do not cross C(n,k).
// Use Dynamic Programming. Use extra space of not more than O(k) number of intals. Do not allocate the whole O(nk) DP table.
// Don't let C(1000,900) take more time than C(1000,500). Time limit may exceed otherwise.
char* intal_bincoeff(unsigned int n, unsigned int k)
{
    if(k > 500)
        k = n-k;
    char *dp[k+1]; //array of k intals
    for(int i = 0;i<=k;++i)
    {
        dp[i] = (char*)malloc((SIZE*10)*sizeof(char));
        strcpy(dp[i], "0");
    }
    char *temp = NULL;
    //Initializing dp array to 0 
    strcpy(dp[0],"1");

    for(int i = 1;i<=n;i++)
    {
        for(int j = (i<k)?i:k ; j > 0; j--) //j starts from the minimum of i or k 
        {
            temp = intal_add(dp[j], dp[j-1]);
            strcpy(dp[j],temp);
            free(temp);
        }
    }
    //free the remaining allocated k memory other than the kth one which holds the final answer
    for(int i = 0;i<k;++i)
    {
        free(dp[i]);
    }
    return dp[k];
}

// Returns the offset of the largest intal in the array.
// Return the smallest offset if there are multiple occurrences.
// 1 <= n <= 1000
int intal_max(char **arr, int n)
{
    int comp_val; //The value obtained from comparison function
    int max = 0;
    for(int i = 1; i< n; ++i)
    {
        comp_val = intal_compare(arr[i], arr[max]);
        if(comp_val == 1)
            max = i;
    }
    return max;
}

// Returns the offset of the smallest intal in the array.
// Return the smallest offset if there are multiple occurrences.
// 1 <= n <= 1000
int intal_min(char **arr, int n)
{
    int comp_val; //The value obtained from comparison function
    int min = 0;
    for(int i = 1; i< n; ++i)
    {
        comp_val = intal_compare(arr[i], arr[min]);
        if(comp_val == -1)
            min = i;
    }
    return min;
}

// Returns the offset of the first occurrence of the key intal in the array.
// Returns -1 if the key is not found.
// 1 <= n <= 1000
int intal_search(char **arr, int n, const char* key)
{
    int found = -1;
    for(int i =0;i<n;++i)
    {
        if(intal_compare(arr[i], key) == 0)
        {
            found = i;
            break;
        }
    }
    return found;
}

// Returns the offset of the first occurrence of the key intal in the SORTED array.
// Returns -1 if the key is not found.
// The array is sorted in nondecreasing order.
// 1 <= n <= 1000
// The implementation should be a O(log n) algorithm.
int intal_binsearch(char **arr, int n, const char* key)
{
    int found = -1;
    int start = 0;
    int end = n-1;
    int m;
    while(end >= start)
    {
        m = (start + end)/2;
        if(intal_compare(arr[m], key) == 0)
        {
            found = m;
            break;
        }
        else if(intal_compare(arr[m],key) == 1)
        {
            end = m-1;
        }
        else{
            start = m+1;
        }
    }
    return found;
}

static int partition(char **arr, int low, int high)
{
    int pivot = low;
    low++;
    while(low <= high)
    {
        while(low <= high && (intal_compare(arr[pivot], arr[low])==1))
            low++;
        while(low <= high && (intal_compare(arr[high], arr[pivot])==1))
            high--;
        if(low < high)
        {
            swap(&arr[low], &arr[high]);
            low++;
            high--;
        }
    } 
    swap(&arr[high], &arr[pivot]);
    return high;  
}

static void quickSort(char **arr, int low, int high)
{
    int pos; //To hold the position of partition element
    if(low < high)
    {
        pos = partition(arr, low, high);

        quickSort(arr, low, pos-1);
        quickSort(arr, pos + 1, high);
    }
}

// Sorts the array of n intals.
// 1 <= n <= 1000
// The implementation should be a O(n log n) algorithm.
void intal_sort(char **arr, int n)
{
    int low = 0;
    int high = n-1;
    quickSort(arr, low, high);
}

// Coin-Row Problem - Dynamic Programming Solution
// There is a row of n coins whose values are some positive integers C[0..n-1].
// The goal is to pick up the maximum amount of money subject to the constraint that
// no two coins adjacent in the initial row can be picked up.
// 1 <= n <= 1000
// The implementation should be O(n) time and O(1) extra space even though the DP table may be of O(n) size.
// Eg: Coins = [10, 2, 4, 6, 3, 9, 5] returns 25
char* coin_row_problem(char **arr, int n)
{
    char *prev = (char*)malloc((SIZE*10)*sizeof(char));
    char *curr = (char*)malloc((SIZE*10)*sizeof(char));
    char *next = (char*)malloc((SIZE*10)*sizeof(char));
    char *temp = (char*)malloc((SIZE*10)*sizeof(char));
    char *temp1 = NULL;
    //Initilaization
    strcpy(prev, "0");
    strcpy(curr, arr[0]);
    for(int i =2; i<=n;++i)
    {
        //if(intal_compare((intal_add(prev, arr[i-1])), curr) == 1) 
        temp1 = intal_add(prev, arr[i-1]);
        strcpy(temp, temp1);
        free(temp1);
        if(intal_compare(temp, curr) == 1)
            strcpy(next, temp);
        else
            strcpy(next, curr);

        strcpy(prev, curr);
        strcpy(curr, next);
    }
    free(prev);
    free(next);
    free(temp);
    return curr;
}
