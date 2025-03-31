// suitable for Arduino NANO board
//DEFINVE X axis
int PULSE_p = 2; //D2
int DIR_p = 3;   //D3
//DEFINVE Y axis
int PULSE1_p = 4; 
int DIR1_p = 5;  
//Trigger pin
int vol = 14;

void setup() {
  pinMode(PULSE_p, OUTPUT);
  pinMode(DIR_p, OUTPUT);
  pinMode(PULSE1_p, OUTPUT);
  pinMode(DIR1_p, OUTPUT);
  Serial.begin(115200);
}

int i = 1600; // for x
int i1 = 0; // for y
int x_dir = 1;
int x1_dir = 1;
int status_on = 0;
unsigned long last_adc_check_time = 0; // 上次检查 ADC 的时间
int value_voltage = 0; // 缓存 ADC 值
bool sendFlag = false; // 标志变量，用于控制是否已经发送过字符

void loop()
{
  //unsigned long current_time = millis(); // 获取当前时间
    //if (current_time - last_adc_check_time > 1000) { // 每隔10ms检查一次
        value_voltage = analogRead(vol); // 更新缓存值,  3.3V is 675;
       // last_adc_check_time = current_time;
    //}
  //Serial.print("The voltage stp is "); Serial.println(value_voltage);

  if (value_voltage > 600 )
    {
      if (!sendFlag) { // 如果尚未发送过字符
        Serial.print("Y"); // 发送字符
        sendFlag = true; // 设置标志变量为 true
      }
      //Serial.print("当前正在正常运行,"); 
      if(i < 3200 && i >= 1600 & x_dir)  { digitalWrite(DIR_p, LOW);   i = i+1;}//Direction n
      else if(i ==  3200 && x_dir)  { digitalWrite(DIR_p, HIGH);   i = i-1; x_dir = 0;}//Direction p
      else if(i > 0 && x_dir == 0)  { digitalWrite(DIR_p, HIGH);   i = i-1; x_dir = 0;}//Direction p
      else if(i == 0 && x_dir == 0) { digitalWrite(DIR_p, LOW);   i = i+1; x_dir = 1;}
      else if(i > 0 &&  i< 1600 && x_dir)  { digitalWrite(DIR_p, LOW);   i = i+1; x_dir = 1;}
      //Serial.print("i = "); Serial.print(i);

      if(i1 < 3200 && x1_dir) {digitalWrite(DIR1_p, LOW); i1 = i1 +1;} //move to motor
      else if(i1 == 3200) {digitalWrite(DIR1_p, HIGH); i1 = i1 -1; x1_dir = 0;}
      else if(x1_dir == 0 && i1> 0) {digitalWrite(DIR1_p, HIGH); i1 = i1 -1; }
      else if(x1_dir == 0 && i1 == 0) {digitalWrite(DIR1_p, LOW); i1 = i1 +1; x1_dir = 1;}
      //Serial.print(", i1 = "); Serial.println(i1);

      // send a pusle, control movement on X direction, 400 pulse/round,  move 8mm in 2 seconds, each pulse period is 2s/(400*8)=2s/3200=2 ms/3.2 = 0.625 ms =625 us
      // send a pusle, control movement on Y direction, 800 pulse/round,  move 4mm in 2 seconds, each pulse period us 2s/(3200)= 625 us
      //digitalWrite(PULSE_p, HIGH); // X motor
      digitalWrite(PULSE1_p, HIGH); // Y motor
      delayMicroseconds(255);//delay 0.2ms
      delayMicroseconds(455);//delay 0.2ms
      
      //digitalWrite(PULSE_p, LOW); 
      digitalWrite(PULSE1_p, LOW); 
      delayMicroseconds(235);//delay 0.2ms
      delayMicroseconds(435);//delay 0.2ms
    }

  if (value_voltage <= 600 ) // operation stop, the chamber come back the prigional position
    {
       //Serial.print("当前正在回归初始位置，"); 
       Serial.print("N");
       if(i< 1600) {
                    digitalWrite(DIR_p, LOW);   
                    i = i+1;
                    // send a pusle, control movement on X direction, 400 pulse/round,  move 8mm in 2 seconds, each pulse period is 2s/(400*8)=2s/3200=2 ms/3.2 = 0.625 ms =625 us
                    digitalWrite(PULSE_p, HIGH);
                    delayMicroseconds(312);//delay 0.2ms
                    digitalWrite(PULSE_p, LOW); 
                    delayMicroseconds(313);//delay 0.2ms
                   }//Direction n
       else if(i > 1600) { 
                    digitalWrite(DIR_p, HIGH);   
                    i = i-1;
                    digitalWrite(PULSE_p, HIGH);
                    delayMicroseconds(312);//delay 0.2ms
                    digitalWrite(PULSE_p, LOW); 
                    delayMicroseconds(313);//delay 0.2ms
                    }//Direction p
        else if(i == 1600) {
                    digitalWrite(PULSE_p, LOW); 
                    delayMicroseconds(313);//delay 0.2ms  
                    }
       //Serial.print("i =");   Serial.print(i); 

       if(i1 != 0) { 
            digitalWrite(DIR_p, HIGH);   
            i1 = i1-1;
            // send a pusle, control movement on Y direction, 800 pulse/round,  move 4mm in 2 seconds, each pulse period us 2s/(3200)= 625 us
            digitalWrite(PULSE1_p, HIGH);
            delayMicroseconds(312);//delay 0.2ms
            digitalWrite(PULSE1_p, LOW); 
            delayMicroseconds(313);//delay 0.2ms
            }//Direction n
       
       else if(i1 == 0)
          { 
            digitalWrite(PULSE1_p, LOW); 
            delayMicroseconds(313);//delay 0.2ms
          }
       
       //Serial.print(", i1 =");   Serial.println(i1); 
    }

}
