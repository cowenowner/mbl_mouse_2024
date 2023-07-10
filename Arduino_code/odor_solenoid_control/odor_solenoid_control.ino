//////////////////////////////////////////////////////
// Code for delivering a sequence of odors to the headfixed mouse
// Cowen and NS&B Students 2023
//////////////////////////////////////////////////////
int n_solenoids = 4;
//byte solenoid_pins[] = {11, 12, 13}; 
byte solenoid_pins[] = {10, 11, 12}; 
byte neutral_odor_pin = 13;

int sequence_length = 6;
byte odor_sequence[] = {1, 2, 3, 1, 2, 3}; // so far, pin/odor port 9 does not work
byte ttl_odor_code_pin = 8;
int ttl_time_ms = 40;  
int port_open_time_ms = 2000;
//int port_closed_time_ms = 50;
int port_to_activate;
int neutral_port_open_time_ms = 3000;

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  Serial.begin(9600);
  Serial.println("STARTED");
  pinMode(neutral_odor_pin, OUTPUT);
  for (int i = 0;i<n_solenoids;i++){
    pinMode(solenoid_pins[i], OUTPUT);
  //  digitalWrite(solenoid_pins[i], LOW);
  }
}

// the loop function runs over and over again forever
void loop() {

    // digitalWrite(ttl_odor_code_pin, HIGH);
    // delay(ttl_time_ms);
    // digitalWrite(ttl_odor_code_pin, LOW);
    // delay(ttl_time_ms);

  for (int i = 0;i<sequence_length;i++){
    port_to_activate = solenoid_pins[odor_sequence[i] - 1];
    Serial.print("Sol ");
    Serial.print(i+1);
    Serial.print(", Pin ");
    Serial.println(port_to_activate);

    // Send a code to Neuropixels to indicate the odor type.
    for (int j = 0;j<odor_sequence[i];j++){
      digitalWrite(ttl_odor_code_pin, HIGH);
      delay(ttl_time_ms);
      digitalWrite(ttl_odor_code_pin, LOW);
      delay(ttl_time_ms);
    }
    // activate the neutral pin
    digitalWrite(port_to_activate, HIGH);
    delay(port_open_time_ms);                     // wait 
    digitalWrite(port_to_activate, LOW);
    //delay(port_closed_time_ms);                   // wait 
     // send out the neutral odor
    digitalWrite(neutral_odor_pin, HIGH);
    delay(neutral_port_open_time_ms);                     // wait 
    digitalWrite(neutral_odor_pin, LOW);
    //delay(port_closed_time_ms);                   // wait 
  }
}
