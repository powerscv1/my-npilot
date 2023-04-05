import copy
import random
import numpy as np
from common.numpy_fast import clip, interp
from cereal import car
from common.realtime import DT_CTRL
from common.conversions import Conversions as CV
from selfdrive.car.hyundai.values import Buttons
from common.params import Params
from selfdrive.controls.lib.drive_helpers import V_CRUISE_MAX, V_CRUISE_MIN, CONTROL_N
from selfdrive.controls.lib.lateral_planner import TRAJECTORY_SIZE
#from selfdrive.controls.lib.longitudinal_mpc_lib.long_mpc import AUTO_TR_CRUISE_GAP
from selfdrive.car.hyundai.values import CAR
from selfdrive.car.isotp_parallel_query import IsoTpParallelQuery
from common.filter_simple import StreamingMovingAverage
from cereal import log


from selfdrive.road_speed_limiter import road_speed_limiter_get_max_speed, road_speed_limiter_get_active, \
  get_road_speed_limiter
#import common.loger as trace1

CREEP_SPEED = 2.3
SYNC_MARGIN = 3.

# do not modify
MIN_SET_SPEED_KPH = V_CRUISE_MIN
MAX_SET_SPEED_KPH = V_CRUISE_MAX

ALIVE_COUNT = [6, 8]
WAIT_COUNT = [12, 13, 14, 15, 16]
AliveIndex = 0
WaitIndex = 0

MIN_CURVE_SPEED = 20. * CV.KPH_TO_MS

EventName = car.CarEvent.EventName

ButtonType = car.CarState.ButtonEvent.Type
ButtonPrev = ButtonType.unknown
ButtonCnt = 0
LongPressed = False

XState = log.LongitudinalPlan.XState

## 국가법령정보센터: 도로설계기준
V_CURVE_LOOKUP_BP = [0., 1./670., 1./560., 1./440., 1./360., 1./265., 1./190., 1./135., 1./85., 1./55., 1./30., 1./15.]
#V_CRUVE_LOOKUP_VALS = [300, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20]
V_CRUVE_LOOKUP_VALS = [300, 120, 110, 100, 90, 80, 70, 60, 50, 45, 35, 30]


class CruiseHelper:

  def __init__(self):
    self.mapValid = 0
    self.v_cruise_kph_apply = 20
    self.v_cruise_kph_backup = 20
    self.curve_speed_last = 255
    self.longActiveUser = 0
    self.preBrakePressed = False
    self.preGasPressedMax = 0.0
    self.gasPressedCount = 0
    self.curvature = 0
    self.position_x = 1000.0
    self.position_y = 300.0
    self.cruiseButtons = 0
    self.userCruisePaused = False
    self.radarAlarmCount = 0
    self.naviSpeedLimitTarget = 30
    self.dRel = 0
    self.vRel = 0
    self.dRelValidCount = 0
    self.trafficState = 0
    self.xState_prev = XState.cruise
    self.xState = XState.cruise
    self.v_ego_kph = 0
    self.v_ego_kph_set = 0
    self.blinker = False
    self.gasPressedFrame = 0
    self.slowSpeedFrameCount = 0
    self.frame = 0
    self.longActiveUserReady = 0
    self.naviSpeed = 255
    self.roadSpeed = 255

    self.active_cam = False
    self.over_speed_limit = False

    self.roadLimitSpeed = 0.0
    self.ndaActive = 0
    self.apilotEventFrame = 0
    self.apilotEventWait = 0
    self.apilotEventPrev = 0

    self.update_params_count = 0
    self.curvatureFilter = StreamingMovingAverage(10)

    self.longCruiseGap = int(Params().get("PrevCruiseGap"))
    self.cruiseSpeedMin = int(Params().get("CruiseSpeedMin"))

    self.autoCurveSpeedCtrlUse = int(Params().get("AutoCurveSpeedCtrlUse"))
    self.autoCurveSpeedIndex = int(Params().get("AutoCurveSpeedIndex"))
    self.autoCurveSpeedFactor = float(int(Params().get("AutoCurveSpeedFactor", encoding="utf8")))*0.01
    self.autoNaviSpeedCtrl = int(Params().get("AutoNaviSpeedCtrl"))
    self.autoNaviSpeedCtrlStart = float(Params().get("AutoNaviSpeedCtrlStart"))
    self.autoNaviSpeedCtrlEnd = float(Params().get("AutoNaviSpeedCtrlEnd"))
    self.autoRoadLimitCtrl = int(Params().get("AutoRoadLimitCtrl", encoding="utf8"))
    self.naviSpeedLimitDecelRate = float(Params().get("NaviSpeedLimitDecelRate", encoding="utf8"))*0.01
    #self.naviDecelMarginDist = float(int(Params().get("NaviDecelMarginDist", encoding="utf8")))
    #self.naviDecelRate = float(int(Params().get("NaviDecelRate", encoding="utf8")))
    self.autoResumeFromGasSpeed = float(int(Params().get("AutoResumeFromGasSpeed", encoding="utf8")))
    self.autoResumeFromGas = Params().get_bool("AutoResumeFromGas")
    self.autoResumeFromBrakeRelease = Params().get_bool("AutoResumeFromBrakeRelease")
    self.autoSyncCruiseSpeedMax = int(Params().get("AutoSyncCruiseSpeedMax"))
    self.autoResumeFromBrakeReleaseDist = float(int(Params().get("AutoResumeFromBrakeReleaseDist", encoding="utf8")))
    self.autoResumeFromBrakeReleaseLeadCar = Params().get_bool("AutoResumeFromBrakeReleaseLeadCar")
    self.autoResumeFromBrakeCarSpeed = float(int(Params().get("AutoResumeFromBrakeCarSpeed", encoding="utf8")))
    self.autoResumeFromBrakeReleaseTrafficSign  = Params().get_bool("AutoResumeFromBrakeReleaseTrafficSign")
    self.longControlActiveSound = int(Params().get("LongControlActiveSound"))
    #self.accelLimitEcoSpeed = float(int(Params().get("AccelLimitEcoSpeed", encoding="utf8")))
    self.autoSpeedUptoRoadSpeedLimit = float(int(Params().get("AutoSpeedUptoRoadSpeedLimit", encoding="utf8"))) / 100.
    #self.accelLimitConfusedModel = int(Params().get("AccelLimitConfusedModel"))
    self.autoSpeedAdjustWithLeadCar = float(int(Params().get("AutoSpeedAdjustWithLeadCar", encoding="utf8"))) / 1.
    self.cruiseButtonMode = int(Params().get("CruiseButtonMode"))
    self.gapButtonMode = int(Params().get("GapButtonMode"))
    self.autoResumeFromGasSpeedMode = int(Params().get("AutoResumeFromGasSpeedMode"))
    self.myDrivingMode = int(Params().get("InitMyDrivingMode"))
    self.mySafeModeFactor = float(int(Params().get("MySafeModeFactor", encoding="utf8"))) / 100. if self.myDrivingMode == 2 else 1.0
    self.liveSteerRatioApply  = float(int(Params().get("LiveSteerRatioApply", encoding="utf8"))) / 100.
    self.autoCancelFromGasMode = int(Params().get("AutoCancelFromGasMode"))
    self.steerActuatorDelay = float(int(Params().get("SteerActuatorDelay", encoding="utf8"))) / 100.
    self.steerActuatorDelayLow = float(int(Params().get("SteerActuatorDelayLow", encoding="utf8"))) / 100.
    self.steerActuatorDelayMid = float(int(Params().get("SteerActuatorDelayMid", encoding="utf8"))) / 100.

  def update_params(self, frame):
    if frame % 20 == 0:
      self.update_params_count += 1
      self.update_params_count = self.update_params_count % 20

      if self.update_params_count == 0:
        self.autoCurveSpeedCtrlUse = int(Params().get("AutoCurveSpeedCtrlUse"))
        self.autoCurveSpeedFactor = float(int(Params().get("AutoCurveSpeedFactor", encoding="utf8")))*0.01
      elif self.update_params_count == 1:
        self.autoNaviSpeedCtrl = int(Params().get("AutoNaviSpeedCtrl"))
        self.autoRoadLimitCtrl = int(Params().get("AutoRoadLimitCtrl", encoding="utf8"))
      elif self.update_params_count == 2:
        self.naviSpeedLimitDecelRate = float(Params().get("NaviSpeedLimitDecelRate", encoding="utf8"))*0.01
        #self.naviDecelMarginDist = float(int(Params().get("NaviDecelMarginDist", encoding="utf8")))
      elif self.update_params_count == 3:
        #self.naviDecelRate = float(int(Params().get("NaviDecelRate", encoding="utf8")))
        self.autoResumeFromGasSpeed = float(int(Params().get("AutoResumeFromGasSpeed", encoding="utf8")))
      elif self.update_params_count == 4:
        self.autoResumeFromGas = Params().get_bool("AutoResumeFromGas")
        self.autoResumeFromBrakeRelease = Params().get_bool("AutoResumeFromBrakeRelease")
      elif self.update_params_count == 5:
        self.autoSyncCruiseSpeedMax = int(Params().get("AutoSyncCruiseSpeedMax"))
        self.autoResumeFromBrakeReleaseDist = float(int(Params().get("AutoResumeFromBrakeReleaseDist", encoding="utf8")))
      elif self.update_params_count == 6:
        self.autoResumeFromBrakeReleaseLeadCar = Params().get_bool("AutoResumeFromBrakeReleaseLeadCar")
        self.autoResumeFromBrakeCarSpeed = float(int(Params().get("AutoResumeFromBrakeCarSpeed", encoding="utf8")))
      elif self.update_params_count == 7:
        self.autoResumeFromBrakeReleaseTrafficSign  = Params().get_bool("AutoResumeFromBrakeReleaseTrafficSign")
        self.longControlActiveSound = int(Params().get("LongControlActiveSound"))
      elif self.update_params_count == 8:
        #self.accelLimitEcoSpeed = float(int(Params().get("AccelLimitEcoSpeed", encoding="utf8")))
        self.autoSpeedUptoRoadSpeedLimit = float(int(Params().get("AutoSpeedUptoRoadSpeedLimit", encoding="utf8"))) / 100.
      elif self.update_params_count == 9:
        #self.accelLimitConfusedModel = int(Params().get("AccelLimitConfusedModel"))
        self.autoSpeedAdjustWithLeadCar = float(int(Params().get("AutoSpeedAdjustWithLeadCar", encoding="utf8"))) / 1.
      elif self.update_params_count == 10:
        self.cruiseButtonMode = int(Params().get("CruiseButtonMode"))
        self.autoResumeFromGasSpeedMode = int(Params().get("AutoResumeFromGasSpeedMode"))
      elif self.update_params_count == 11:
        #self.myDrivingMode = int(Params().get("InitMyDrivingMode")) #초기에 한번만 읽어옴...
        self.mySafeModeFactor = float(int(Params().get("MySafeModeFactor", encoding="utf8"))) / 100. if self.myDrivingMode == 2 else 1.0
        self.liveSteerRatioApply  = float(int(Params().get("LiveSteerRatioApply", encoding="utf8"))) / 100.
      elif self.update_params_count == 12:
        self.autoCancelFromGasMode = int(Params().get("AutoCancelFromGasMode"))
        self.gapButtonMode = int(Params().get("GapButtonMode"))
      elif self.update_params_count == 13:
        self.steerActuatorDelay = float(int(Params().get("SteerActuatorDelay", encoding="utf8"))) / 100.
        self.steerActuatorDelayLow = float(int(Params().get("SteerActuatorDelayLow", encoding="utf8"))) / 100.
        self.steerActuatorDelayMid = float(int(Params().get("SteerActuatorDelayMid", encoding="utf8"))) / 100.
      elif self.update_params_count == 14:
        self.cruiseSpeedMin = int(Params().get("CruiseSpeedMin"))
      elif self.update_params_count == 15:
        self.autoNaviSpeedCtrlStart = float(Params().get("AutoNaviSpeedCtrlStart"))
        self.autoNaviSpeedCtrlEnd = float(Params().get("AutoNaviSpeedCtrlEnd"))
      elif self.update_params_count == 16:
        self.autoCurveSpeedIndex = int(Params().get("AutoCurveSpeedIndex"))


  def getSteerActuatorDelay(self, v_ego):
    v_ego_kph = v_ego * 3.6

    return interp(v_ego_kph, [0, 50, 200], [self.steerActuatorDelayLow, self.steerActuatorDelayMid, self.steerActuatorDelay])

  @staticmethod
  def get_lead(sm):

    radar = sm['radarState']
    if radar.leadOne.status and radar.leadOne.radar:
      return radar.leadOne

    return None

# ajouatom
# navi 운영방법
# - 속도제한신호(mapValid)가 들어오면 현재크루즈속도(v_cruise_kph) 백업(v_cruise_kph_backup)
# - 속도제한신호가 나가면 백업된 크르즈속도 복원
# - 운행중 사용자가 속도변경을 하면?
#   - 버튼을 눌렀을때: resume버튼시 백업속도로 변경, 현재크루즈속도를 백업하도록 함.
#   - 엑셀을 밟았을때: 현재크루즈속도변경, 현재크르주속도를 백업하도록 함.
#   - 브레이크를 밟았을때: 현재크루즈속도 변경, resume시 백업속도로 주행
# - 속도제한신호가 들어오면
# - 구간단속은???? 몰라..
# - 크루즈 운용방법(HW)
#   - 브레이크를 밟으면: Cruise disable... 
#   - 브레이클를 밟았다 떼면: Cruise disable... 유지
#   - 엑셀을 밟으면: 30Km/h이상되면 auto Resume
#   - 엑셀을 밟았다 떼면: zz

  def cruise_control(self, controls, CS, active_mode=0):  #active_mode => -3(OFF auto), -2(OFF brake), -1(OFF user), 0(OFF), 1(ON user), 2(ON gas), 3(ON auto)
    if controls.enabled:
      if active_mode > 0 and controls.CC.longEnabled:
        if self.longActiveUser <= 0:
          controls.LoC.reset(v_pid=CS.vEgo)
        if self.longControlActiveSound >= 2 and self.longActiveUser != active_mode:
          controls.events.add(EventName.cruiseResume)
        self.longActiveUser = active_mode
        self.userCruisePaused = False
        if self.gapButtonMode == 3:
          self.longCruiseGap = 4
      elif active_mode <= 0:
        if self.longActiveUser != active_mode and self.longControlActiveSound >= 2:
          #controls.events.add(EventName.cruisePaused)
          pass
        self.longActiveUser = active_mode

  @staticmethod
  def pcmCruiseControl(v_cruise_kph_last, v_cruise_kph):
    pass

  def get_lead_rel(self, controls):    
    lead = self.get_lead(controls.sm)
    dRel = lead.dRel if lead is not None else 0
    vRel = lead.vRel if lead is not None else 0
    return dRel, vRel

  def update_cruise_buttons(self, enabled, buttonEvents, v_cruise_kph, metric):
    global ButtonCnt, LongPressed, ButtonPrev

    button_speed_up_diff = 1
    button_speed_dn_diff = 10 if self.cruiseButtonMode in [3, 4] else 1

    button_type = 0
    if enabled:
      if ButtonCnt:
        ButtonCnt += 1
      for b in buttonEvents:
        if b.pressed and not ButtonCnt and (b.type == ButtonType.accelCruise or b.type == ButtonType.decelCruise or b.type == ButtonType.gapAdjustCruise):
          ButtonCnt = 1
          ButtonPrev = b.type
        elif not b.pressed and ButtonCnt:
          if not LongPressed and b.type == ButtonType.accelCruise:
            v_cruise_kph += button_speed_up_diff if metric else button_speed_up_diff * CV.MPH_TO_KPH
            button_type = ButtonType.accelCruise
          elif not LongPressed and b.type == ButtonType.decelCruise:
            v_cruise_kph -= button_speed_dn_diff if metric else button_speed_dn_diff * CV.MPH_TO_KPH
            button_type = ButtonType.decelCruise
          elif not LongPressed and b.type == ButtonType.gapAdjustCruise:
            if self.gapButtonMode == 0:
              self.longCruiseGap = 1 if self.longCruiseGap == 4 else self.longCruiseGap + 1
            elif self.gapButtonMode == 1:
              self.longCruiseGap = 1 if self.longCruiseGap == 5 else self.longCruiseGap + 1
            elif self.gapButtonMode == 2:
              self.longCruiseGap = 4 if self.longCruiseGap == 5 else self.longCruiseGap + 1
            elif self.gapButtonMode == 3:
              self.longCruiseGap = 5

            button_type = ButtonType.gapAdjustCruise
            #Params().put("PrevCruiseGap", str(self.longCruiseGap))

          LongPressed = False
          ButtonCnt = 0
      if ButtonCnt > 40:
        LongPressed = True
        V_CRUISE_DELTA = 10
        if ButtonPrev == ButtonType.accelCruise:
          v_cruise_kph += V_CRUISE_DELTA - v_cruise_kph % V_CRUISE_DELTA
          button_type = ButtonType.accelCruise
          ButtonCnt %= 40
        elif ButtonPrev == ButtonType.decelCruise:
          v_cruise_kph -= V_CRUISE_DELTA - -v_cruise_kph % V_CRUISE_DELTA
          button_type = ButtonType.decelCruise
          ButtonCnt %= 40
        elif ButtonPrev == ButtonType.gapAdjustCruise:
          button_type = ButtonType.gapAdjustCruise
          ButtonCnt = 0
    v_cruise_kph = clip(v_cruise_kph, self.cruiseSpeedMin, MAX_SET_SPEED_KPH)
    return button_type, LongPressed, v_cruise_kph

  def update_speed_navi(self, CS, controls, v_cruise_kph):
    navi = CS.naviSafetyInfo
    #controls.debugText1 = 'S{}/{},D{}/{},N{}'.format(navi.sign, navi.speed2, navi.dist1, navi.dist2, navi.speedLimit)
    v_ego_kph = int(CS.vEgo * CV.MS_TO_KPH + 0.5) + 3.0
    v_kph_apply = v_cruise_kph
    if self.naviSpeedLimitDecelRate > 0.0 and navi.speedLimit > 0 and navi.speedLimit < 255 and v_kph_apply > navi.speedLimit:
      v_kph_apply = min(v_ego_kph, self.naviSpeedLimitTarget)
      self.naviSpeedLimitTarget = max(navi.speedLimit, v_kph_apply - self.naviSpeedLimitDecelRate * 2.5 * DT_CTRL)
    else:
      self.naviSpeedLimitTarget = v_cruise_kph
      return 0
    return self.naviSpeedLimitTarget

  def update_speed_nda(self, CS, controls):
    clu11_speed = CS.vEgoCluster * CV.MS_TO_KPH
    road_speed_limiter = get_road_speed_limiter()
    self.ndaActive = 1 if road_speed_limiter_get_active() > 0 else 0
    apply_limit_speed, road_limit_speed, left_dist, first_started, max_speed_log = \
      road_speed_limiter.get_max_speed(clu11_speed, True, self.autoNaviSpeedCtrlStart, self.autoNaviSpeedCtrlEnd) #self.is_metric)

    self.active_cam = road_limit_speed > 0 and left_dist > 0

    if road_speed_limiter.roadLimitSpeed is not None:
      camSpeedFactor = clip(road_speed_limiter.roadLimitSpeed.camSpeedFactor, 1.0, 1.1)
      self.over_speed_limit = road_speed_limiter.roadLimitSpeed.camLimitSpeedLeftDist > 0 and \
                              0 < road_limit_speed * camSpeedFactor < clu11_speed + 2
    else:
      self.over_speed_limit = False

    #str1 = 'applyLimit={},speedLimit={},leftDist={}'.format(apply_limit_speed, road_limit_speed, left_dist)
    #controls.debugText1 = str1
    self.roadLimitSpeed = controls.sm['roadLimitSpeed'].roadLimitSpeed

    return clip(apply_limit_speed, 0, MAX_SET_SPEED_KPH), clip(self.roadLimitSpeed, 30, MAX_SET_SPEED_KPH)

  def apilot_curve(self, CS, controls):
    curvatures = controls.sm['lateralPlan'].curvatures
    turnSpeed = 300
    if len(curvatures) == CONTROL_N:
      curvature = abs(self.curvatureFilter.process(curvatures[self.autoCurveSpeedIndex]))  * self.autoCurveSpeedFactor
      if abs(curvature) > 0.001:
        turnSpeed = interp(curvature, V_CURVE_LOOKUP_BP, V_CRUVE_LOOKUP_VALS)
        turnSpeed = clip(turnSpeed, MIN_CURVE_SPEED, MAX_SET_SPEED_KPH)
      else:
        turnSpeed = 300
    else:
      self.curvatureFilter.set(0.0)

    controls.debugText1 = 'CURVE={:5.1f},curvature={:5.4f}'.format(turnSpeed, curvature)
    self.turnSpeed_prev = turnSpeed
    return turnSpeed

  def v_cruise_speed_up(self, v_cruise_kph, roadSpeed):
    if v_cruise_kph < roadSpeed:
      v_cruise_kph = roadSpeed
    else:
      for speed in range (40, 150, 10):
        if v_cruise_kph < speed:
          v_cruise_kph = speed
          break
    return clip(v_cruise_kph, self.cruiseSpeedMin, MAX_SET_SPEED_KPH)

  def send_apilot_event(self, controls, eventName, waiting = 20):
    #if eventName != self.apilotEventPrev or (controls.sm.frame - self.apilotEventFrame)*DT_CTRL > self.apilotEventWait: 
    if (controls.sm.frame - self.apilotEventFrame)*DT_CTRL > self.apilotEventWait: # 시끄러..
       controls.events.add(eventName)
       self.apilotEventFrame = controls.sm.frame
       self.apilotEventPrev = eventName
       self.apilotEventWait = waiting

  ######################
  # Gas(엑셀) 놓는 시점
    #### CruiseOFF상태이면
      #  1. 가속페달 CruiseON (autoResumeFromGas) & 신호적색아님 & (autoResumeFromGasSpeed보다 빠거나 60%이상 밟으면
    #### CruiseON상태이면
      #  1. softHold상태: cruiseOFF: 엑셀로 밟으면 크루즈해제
      #  2. 신호감지감속중: cruiseOFF: 신호감지감속이 맘에 안드는 상태, 가속페달을 밟으면 해제
      #  3. 저속주행: cruiseOFF(autoResumeFromGasSpeed 이하): 조건(autoCancelFromGasMode)에 따라 선행차의 유무에 따라 크루즈 해제
      #  4. 설정속도보다 느림 & 커브 & 선행차30M이내: 크루즈속도 현재속도셋
      #  5. 페달을 0.6초이내 뗀경우: 속도증가 autoSyncCruiseSpeedMax까지: 가속페달로 속도를 증가시킴
      #  6. 크루즈속도보다 높을때: 크루즈속도 현재속도셋 : autoSyncCruiseSpeedMax까지
  #####################
  def check_gas_cruise_on(self, CS, v_cruise_kph):
    resume_cond = abs(CS.steeringAngleDeg) < 20 # and not CS.steeringPressed
    longActiveUser = self.longActiveUser
    v_cruise_kph_backup = self.v_cruise_kph_backup
    
    #### CruiseOFF상태이면
    if self.longActiveUser <= 0:
      #  1. 가속페달 CruiseON (autoResumeFromGas) & 신호적색아님 & (autoResumeFromGasSpeed보다 빠거나 60%이상 밟으면
      #    - autoResumeFromGasSpeedMode에 따라 속도 설정(기존속도, 현재속도)
      if ((resume_cond and (self.v_ego_kph >= self.autoResumeFromGasSpeed)) or CS.gas >= 0.6) and self.trafficState != 1 and self.autoResumeFromGas:
        if self.autoResumeFromGasSpeedMode == 0: #현재속도로 세트
          # 60% 이상 밟으면...
          if self.preGasPressedMax >= 0.6:
            v_cruise_kph = self.v_cruise_kph_backup # 60%이상 GAS를 밟으면..
          else:
            v_cruise_kph = self.v_ego_kph_set  # 현재속도로 세트~
        elif self.autoResumeFromGasSpeedMode == 1:   #기존속도
            v_cruise_kph = self.v_cruise_kph_backup 
        elif self.autoResumeFromGasSpeedMode == 2:   #레이더가 검출될때만 기존속도..
          if self.dRel > 0:
            v_cruise_kph = self.v_cruise_kph_backup 
          else:
            v_cruise_kph = self.v_ego_kph_set  # 현재속도로 세트~
        if self.longCruiseGap != 5: 
          longActiveUser = 3
    #### CruiseON상태이면
    else:
      #  1. softHold상태: cruiseOFF: 엑셀로 밟으면 크루즈해제
      if self.xState == XState.softHold:
        longActiveUser = -2
      #  2. 신호감지감속중: cruiseOFF: 신호감지감속이 맘에 안드는 상태, 가속페달을 밟으면 해제
      elif self.xState in [XState.e2eStop, XState.e2eCruise] and self.v_ego_kph < v_cruise_kph and self.trafficState == 1: #controls.v_future*CV.MS_TO_KPH < v_ego_kph * 0.6: 
        longActiveUser = -2
      #  3. 저속주행: cruiseOFF(autoResumeFromGasSpeed 이하): 조건(autoCancelFromGasMode)에 따라 선행차의 유무에 따라 크루즈 해제
      elif self.v_ego_kph < self.autoResumeFromGasSpeed:
        if self.autoCancelFromGasMode == 1:
          longActiveUser = -2
        elif self.autoCancelFromGasMode == 2 and self.dRel==0: # mode:2일때는 선행차가 없을때만
          longActiveUser = -2
      #  4. 설정속도보다 느림 & 커브 & 선행차30M이내: 크루즈속도 현재속도셋
      elif self.v_ego_kph < v_cruise_kph and abs(CS.steeringAngleDeg) > 7.0 and 0 < self.dRel < 30: 
          v_cruise_kph = self.v_ego_kph_set
      #  5. 페달을 0.6초이내 뗀경우: 속도증가 autoSyncCruiseSpeedMax까지: 가속페달로 속도를 증가시킴
      elif self.gasPressedCount * DT_CTRL < 0.6:
        if self.v_ego_kph > 30.0 and v_cruise_kph < self.autoSyncCruiseSpeedMax:  
          v_cruise_kph = self.v_cruise_speed_up(v_cruise_kph, self.roadSpeed)
          if self.autoSyncCruiseSpeedMax > 0 and v_cruise_kph > self.autoSyncCruiseSpeedMax:
            v_cruise_kph = self.autoSyncCruiseSpeedMax
          self.v_cruise_kph_backup = v_cruise_kph
      #  6. 크루즈속도보다 높을때: 크루즈속도 현재속도셋 : autoSyncCruiseSpeedMax까지
      elif self.v_ego_kph > v_cruise_kph and self.autoSyncCruiseSpeedMax > self.autoResumeFromGasSpeed:
        if self.autoResumeFromGasSpeed < self.v_ego_kph < self.autoSyncCruiseSpeedMax: # 오토크루즈 ON속도보다 높고, 130키로보다 작을때만 싱크
          v_cruise_kph = self.v_ego_kph_set
          self.v_cruise_kph_backup = v_cruise_kph #가스로 할땐 백업

    return longActiveUser, v_cruise_kph, v_cruise_kph_backup

  #########################
  # 브레이크해제시: 크루즈ON
  # 설정: 브레이크해제 크루즈ON사용(GAP 5일때는 안함)
  # 설정: 브레이크해제 크루즈ON:주행중,선행차: 거리
  # 설정: 브레이크해제 크루즈ON:정지,선행차
  # 설정: 브레이크해제 크루즈ON:주행중,속도: 속도
  # 설정: 브레이크해제 크루즈ON:주행중,신호
  # 1. softHold
  # 2. 정지할때 : 20km/h이하 속도
  #   - 주차장에서 정지출발을 자주한다. -> 저속주행10초이상 주행시 가속페달후 5초 이내에 브레이크밟음 -> cruiseOFF 유지
  #   - 끼워들기,좌/우회전 하려고 정지하였다. 20M이내 차량검출됨.  깜박이는 켜저있음.. -> cruiseOFF 유지
  #   - 앞의 차량을 발견하고 정지하였다. 10M -> cruiseON
  #   - 신호등앞에서 깜박이를 켜고 정지하였다. -> 신호오류(Pause) 대기(버튼출발) -> cruiseON (일단 OFF, 추후코딩할것!)
  #   - 신호등앞에서 깜박이 없이 정지하였다. -> 직진 신호대기(신호감지) -> cruiseON
  # 3. 속도를 줄이고 싶을때 
  #   - 멀리 차량이 있어 속도를 줄였다. -> 설정거리(20)이상에서(가까울땐 안함)  -> cruiseON
  #       (왜 가까울땐 안하는지 기억이 안남)
  #   - 멀리 신호가 적색이어서 속도를 줄였다. -> 70키로 이하에서 적색이면 -> cruiseON
  #   - 너무빠르다고 생각하여 속도를 줄였다.  -> 설정속도(40) 이상에서 -> cruiseON
  #
  #
  def check_brake_cruise_on(self, CS, v_cruise_kph):
    longActiveUser = self.longActiveUser
    resume_cond = abs(CS.steeringAngleDeg) < 20 # and not CS.steeringPressed
    v_cruise_kph_backup = self.v_cruise_kph_backup

    # 정지상태, 소프트홀드일때 크루즈 ON
    if self.v_ego_kph < 5.0 and self.xState == XState.softHold and self.longCruiseGap != 5:
      longActiveUser = 3
    # 브레이크해제 켜지고, 크루즈갭이 5가 아닌경우에만 작동.
    elif self.autoResumeFromBrakeRelease and self.longCruiseGap != 5: # 브레이크 해제에 대한 크루즈 ON
      gasTime = (self.frame - self.gasPressedFrame)*DT_CTRL
      # 저속 정지.
      if self.v_ego_kph < 20.0:
        gasWaitTime = 5.0 if self.slowSpeedFrameCount*DT_CTRL > 10.0 else 0.0
        if resume_cond and gasTime >= gasWaitTime:  # 5초이상 이전에 gas를 밟았을때만..(저속주행시에만)
          if 0 < self.dRel < 20.0 and self.blinker: # 앞에 차가 있고 깜박이 켜져있음.
            pass
          elif 0 < self.dRel < 10.0: #앞에 차가있음.
            if self.autoResumeFromBrakeReleaseLeadCar:
              longActiveUser = 3
          elif self.dRel == 0 and self.trafficState==1: #선행차없이 신호검출정지
            if self.blinker: # 깜박이켜면 일단 패스
              pass
            else:  # 깜박이가 꺼져있으면 직진이라 가정
              if self.autoResumeFromBrakeReleaseTrafficSign:
                longActiveUser = 3  
      # 주행
      else:
        if 0 < self.dRel:   # 전방에 차량이 있는경우
          if self.dRel > self.autoResumeFromBrakeReleaseDist: ## 설정값 이상의 거리에서만 작동함... 가까울때는 왜 안하게 했지? 
            longActiveUser = 3
        elif self.trafficState == 1:  # 신호감지된경우
          if self.v_ego_kph < 70.0 and self.autoResumeFromBrakeReleaseTrafficSign:  #속도가 70키로 미만이면 
            longActiveUser = 3
          else:        #속도가 빠르면... pass
            pass
        else:   #그냥 감속한경우, 현재속도세트
          if self.v_ego_kph >= self.autoResumeFromBrakeCarSpeed and self.autoResumeFromBrakeCarSpeed > 0:
            longActiveUser = 3
            v_cruise_kph = self.v_ego_kph_set

    return longActiveUser, v_cruise_kph, v_cruise_kph_backup

  # 1. 소프트홀드
  # 2. 정지, 전방차량이 5km/h, 10M이내의 경우 : 선행차뒤에 정지하고 안전한 거리에서 브레이크를 떼면...
  # 3. 전방차량X, 속도가 40키로 이상
  # 4. 직진, 70키로 미만, 신호정지신호, 깜박이 OFF
  # 5. 직진, 전방차량 20M(설정)이상, 깜박이 OFF
  #########################
  def check_brake_cruise_on_old(self, CS, v_cruise_kph):
    longActiveUser = self.longActiveUser
    resume_cond = abs(CS.steeringAngleDeg) < 20 # and not CS.steeringPressed
    v_cruise_kph_backup = self.v_cruise_kph_backup

    # 정지상태, 소프트홀드일때 크루즈 ON
    if self.v_ego_kph < 5.0 and self.xState == XState.softHold and self.longCruiseGap != 5:
      longActiveUser = 3
    # 브레이크해제 켜지고, 크루즈갭이 5가 아닌경우에만 작동.
    elif self.autoResumeFromBrakeRelease and self.longCruiseGap != 5: # 브레이크 해제에 대한 크루즈 ON
      # 주행중, 전방차량이 20M(변수)이상이면
      if resume_cond and self.v_ego_kph >= 3.0 and 0 < self.autoResumeFromBrakeReleaseDist < self.dRel and not self.blinker: #CS.rightBlinker == False:
        v_cruise_kph = self.v_ego_kph_set  # 현재속도로 세트~
        longActiveUser = 3
      # 핸들5도미만, 70km/h미만, 신호정지신호, 깜박이 OFF
      elif abs(CS.steeringAngleDeg) < 5.0 and self.v_ego_kph < 70.0 and self.trafficState == 1 and self.autoResumeFromBrakeReleaseTrafficSign and not self.blinker:
        v_cruise_kph = self.v_ego_kph_set  
        longActiveUser = 3
      # 전방차량이 없고, 속도가 40km/h(변수) 이상인경우
      elif self.dRel==0 and self.v_ego_kph >= self.autoResumeFromBrakeCarSpeed and self.autoResumeFromBrakeCarSpeed > 0:
        v_cruise_kph = self.v_ego_kph_set
        longActiveUser = 3
      # 정지중, 전방차량이 10M이내인경우
      elif self.v_ego_kph < 5.0 and 2 < self.dRel < 10 and self.autoResumeFromBrakeReleaseLeadCar:
        longActiveUser = 3
    return longActiveUser, v_cruise_kph, v_cruise_kph_backup

  def button_control(self, enabled, controls, CS, v_cruise_kph, buttonEvents, metric):
    longActiveUser = self.longActiveUser
    v_cruise_kph_backup = self.v_cruise_kph_backup
    button,buttonLong,buttonSpeed = self.update_cruise_buttons(enabled, buttonEvents, v_cruise_kph, metric)

    ##### Cruise Button 처리...
    if buttonLong:
      if button in [ButtonType.accelCruise, ButtonType.decelCruise]:
        v_cruise_kph = buttonSpeed
        self.v_cruise_kph_backup = v_cruise_kph #버튼으로할땐 백업
      elif button == ButtonType.gapAdjustCruise:  ##안먹네.... 나중에 보자~
        #myDrivingMode = int(Params().get("MyDrivingMode"))
        self.myDrivingMode = self.myDrivingMode + 1 if self.myDrivingMode < 4 else 1
        #Params().put("MyDrivingMode", str(myDrivingMode))
    else:
      self.cruiseButtons = button
      if button == ButtonType.accelCruise:   
        controls.cruiseButtonCounter += 1
        if self.longActiveUser <= 0:
          longActiveUser = 1
          v_cruise_kph = max(v_cruise_kph, self.v_cruise_kph_backup, self.v_ego_kph_set) #브레이크를 밟기전 속도로 복원..
        else:
          if self.xState == XState.softHold:
            longActiveUser = 1
          #elif self.xState == XState.e2eStop:
          #  pass
          elif self.cruiseButtonMode in [1,2]:
            v_cruise_kph = self.v_cruise_speed_up(v_cruise_kph, self.roadSpeed)
          else:
            v_cruise_kph = buttonSpeed
        self.v_cruise_kph_backup = v_cruise_kph #버튼으로할땐 백업
      elif button == ButtonType.decelCruise:
        controls.cruiseButtonCounter -= 1
        if self.longActiveUser <= 0:
          v_cruise_kph = self.v_ego_kph_set  ## 현재속도도 크루즈세트
          longActiveUser = 1
        else:
          if self.xState == XState.softHold:
            longActiveUser = 1
          if CS.gasPressed and v_cruise_kph < self.v_ego_kph_set:
            v_cruise_kph = self.v_ego_kph_set
            self.v_cruise_kph_backup = v_cruise_kph #버튼으로할땐 백업
          elif self.xState == XState.softHold:
            pass
          elif self.xState == XState.e2eStop and self.v_ego_kph < 5: #5km/h 미만, 신호감속중.. (-)를 누르면 크루즈해제... 이러면 설설가겠지? 다시누르면 정지..
            v_cruise_kph = 3
            longActiveUser = -1
            pass
          elif v_cruise_kph > self.v_ego_kph_set+2 and self.cruiseButtonMode in [1,2]:
            v_cruise_kph = self.v_ego_kph_set
            self.v_cruise_kph_backup = v_cruise_kph #버튼으로할땐 백업
          else:
            if self.cruiseButtonMode==2:
              self.userCruisePaused = True
              longActiveUser = -1
            else:
              v_cruise_kph = buttonSpeed
              self.v_cruise_kph_backup = v_cruise_kph
      elif button == ButtonType.gapAdjustCruise and self.gapButtonMode == 3:
        if self.longActiveUser > 0: # and self.gapButtonMode == 3:
          longActiveUser = -1
        else:
          longActiveUser = 1
          v_cruise_kph = max(v_cruise_kph, self.v_cruise_kph_backup, self.v_ego_kph_set) #브레이크를 밟기전 속도로 복원..

    return longActiveUser, v_cruise_kph, v_cruise_kph_backup

  def update_v_cruise_apilot(self, v_cruise_kph, buttonEvents, enabled, metric, controls, CS):
    longActiveUser = self.longActiveUser
    self.frame = controls.sm.frame
    self.update_params(self.frame)
    self.naviSpeed, self.roadSpeed = self.update_speed_nda(CS, controls)
    carNaviSpeed = self.update_speed_navi(CS, controls, v_cruise_kph)
    
    curveSpeed = 255
    if self.autoCurveSpeedCtrlUse > 0:
      curveSpeed = self.apilot_curve(CS, controls)

    self.v_ego_kph = int(CS.vEgo * CV.MS_TO_KPH + 0.5) + 2.0 #실제속도가 v_cruise_kph보다 조금 빨라 2을 더함.
    self.v_ego_kph_set = clip(self.v_ego_kph, self.cruiseSpeedMin, MAX_SET_SPEED_KPH)
    self.xState = controls.sm['longitudinalPlan'].xState
    dRel, vRel = self.get_lead_rel(controls)
    leadCarSpeed = self.v_ego_kph + vRel*CV.MS_TO_KPH
    
    #if 2 < dRel < 20:
    #  self.dRelValidCount += 1
    #elif dRel != 0 or self.dRelValidCount < 50 or abs(CS.steeringAngleDeg)<30:
    #  self.dRelValidCount = 0
    #else:
    #  self.dRelValidCount = 0
    #  if self.radarAlarmCount == 0:
    #    self.radarAlarmCount = 2000
    #    #v_cruise_kph = min(self.v_ego_kph_set, v_cruise_kph) # 레이더가 갑자기 사라지는 경우 현재속도로 세트함.
    
    trafficState = (controls.sm['longitudinalPlan'].trafficState % 100)
    trafficError = controls.sm['longitudinalPlan'].trafficState >= 1000
    if self.longActiveUser>0:
      if self.xState != self.xState_prev and self.xState == XState.softHold:
        controls.events.add(EventName.autoHold)
      if self.xState == XState.softHold and self.trafficState != 2 and trafficState == 2:
        self.send_apilot_event(controls, EventName.trafficSignChanged)
        #self.radarAlarmCount = 2000 if self.radarAlarmCount == 0 else self.radarAlarmCount
      elif self.xState == XState.e2eCruise and self.trafficState != 2 and trafficState == 2 and CS.vEgo < 0.1:
        controls.events.add(EventName.trafficSignGreen)
      elif self.xState == XState.e2eStop and self.xState_prev in [XState.e2eCruise, XState.lead]: # and self.longControlActiveSound >= 2:
        self.send_apilot_event(controls, EventName.trafficStopping, 20.0)
      elif trafficError:
        self.send_apilot_event(controls, EventName.trafficError, 10.0)


    self.trafficState = trafficState
    self.dRel = dRel
    self.vRel = vRel
    self.radarAlarmCount = self.radarAlarmCount - 1 if self.radarAlarmCount > 0 else 0

    self.blinker = CS.rightBlinker or CS.leftBlinker

    longActiveUser, v_cruise_kph, self.v_cruise_kph_backup = self.button_control(enabled, controls, CS, v_cruise_kph, buttonEvents, metric)
    if controls.enabled:      

      if CS.brakePressed:
        longActiveUser = -2
        if not self.preBrakePressed:
          self.v_cruise_kph_backup = v_cruise_kph
        self.longActiveUserReady,temp,temp = self.check_brake_cruise_on(CS, v_cruise_kph)
      elif CS.gasPressed:  
        self.longActiveUserReady,temp,temp = self.check_gas_cruise_on(CS, v_cruise_kph)
      elif not CS.gasPressed and self.gasPressedCount > 2:
        longActiveUser,v_cruise_kph,self.v_cruise_kph_backup = self.check_gas_cruise_on(CS, v_cruise_kph)
      elif not CS.brakePressed and self.preBrakePressed:
        longActiveUser,v_cruise_kph,self.v_cruise_kph_backup = self.check_brake_cruise_on(CS, v_cruise_kph)
      elif self.userCruisePaused:
        if self.v_ego_kph > 3.0 and self.dRel > 0 and self.vRel < 0:          
          v_cruise_kph = self.v_ego_kph_set
          longActiveUser = 3
        elif self.v_ego_kph > 20.0 and self.xState == XState.e2eStop and abs(self.position_y) < 3.0:
          v_cruise_kph = self.v_ego_kph_set
          longActiveUser = 3
        pass

      self.cruise_control(controls, CS, longActiveUser)

      ###### 크루즈 속도제어~~~
      self.v_cruise_kph_apply = v_cruise_kph

      ###### leadCar 관련 속도처리
      roadSpeed1 = self.roadSpeed * self.autoSpeedUptoRoadSpeedLimit
      if v_cruise_kph < roadSpeed1 and self.dRel > 0 and self.vRel > 0 and self.autoSpeedUptoRoadSpeedLimit > 0:
        if leadCarSpeed > v_cruise_kph:
          v_cruise_kph = max(v_cruise_kph, min(leadCarSpeed, roadSpeed1))
          self.v_cruise_kph_apply = v_cruise_kph
      elif self.autoSpeedAdjustWithLeadCar > 0.0 and self.dRel > 0:
        leadCarSpeed1 = max(leadCarSpeed + self.autoSpeedAdjustWithLeadCar, 30)
        if leadCarSpeed1 < v_cruise_kph:
          self.v_cruise_kph_apply = leadCarSpeed1
      #controls.debugText1 = 'LC={:3.1f},{:3.1f},RS={:3.1f},SS={:3.1f}'.format( leadCarSpeed, vRel*CV.MS_TO_KPH, self.roadSpeed, self.v_cruise_kph_apply)      

      ###### naviSpeed, roadSpeed, curveSpeed처리
      if self.autoNaviSpeedCtrl > 0 and self.naviSpeed > 0:
        if self.naviSpeed < self.v_cruise_kph_apply:
          #self.send_apilot_event(controls, EventName.speedDown, 60.0)  #시끄러..
          pass
        self.v_cruise_kph_apply = min(self.v_cruise_kph_apply, self.naviSpeed)
        self.ndaActive = 2 if self.ndaActive == 1 else 0
      elif self.autoNaviSpeedCtrl > 1 and carNaviSpeed > 0:
        self.v_cruise_kph_apply = min(self.v_cruise_kph_apply, carNaviSpeed)
        self.ndaActive = 2
      if self.roadSpeed > 30 and False: # 로드스피드리밋 사용안함..
        if self.autoRoadLimitCtrl == 1:
          self.v_cruise_kph_apply = min(self.v_cruise_kph_apply, self.roadSpeed)
        elif self.autoRoadLimitCtrl == 2:
          self.v_cruise_kph_apply = min(self.v_cruise_kph_apply, self.roadSpeed)
      if self.autoCurveSpeedCtrlUse > 0:
        if curveSpeed < self.v_cruise_kph_apply and self.longActiveUser > 0:
          #self.send_apilot_event(controls, EventName.speedDown, 60.0)
          pass
        self.v_cruise_kph_apply = min(self.v_cruise_kph_apply, curveSpeed)
    else: #not enabled
      self.v_cruise_kph_backup = v_cruise_kph #not enabled

    self.preBrakePressed = CS.brakePressed
    self.xState_prev = self.xState
    if self.v_ego_kph < 20.0:
      self.slowSpeedFrameCount += 1
    else:
      self.slowSpeedFrameCount = 0

    if CS.gasPressed:
      self.gasPressedFrame = self.frame
      self.gasPressedCount += 1
      if CS.gas > self.preGasPressedMax:
        self.preGasPressedMax = CS.gas
      #controls.debugText1 = 'GAS: {:3.1f}/{:3.1f}={:3.1f}'.format(CS.gas*100., self.preGasPressedMax*100., self.gasPressedCount * DT_CTRL)
    else:
      self.preGasPressedMax = 0.0
      self.gasPressedCount = 0
    return v_cruise_kph

def enable_radar_tracks(CP, logcan, sendcan):
  # START: Try to enable radar tracks
  print("Try to enable radar tracks")  
  # if self.CP.openpilotLongitudinalControl and self.CP.carFingerprint in [HYUNDAI_CAR.SANTA_FE_2022]:
  if CP.openpilotLongitudinalControl: # and CP.carFingerprint in [CAR.SANTA_FE, CAR.SANTA_FE_HEV_2022, CAR.NEXO]:
    rdr_fw = None
    rdr_fw_address = 0x7d0 #일부차량은 다름..
    if True:
      for i in range(10):
        print("O yes")
      try:
        for i in range(40):
          try:
            query = IsoTpParallelQuery(sendcan, logcan, CP.sccBus, [rdr_fw_address], [b'\x10\x07'], [b'\x50\x07'], debug=True)
            for addr, dat in query.get_data(0.1).items(): # pylint: disable=unused-variable
              print("ecu write data by id ...")
              new_config = b"\x00\x00\x00\x01\x00\x01"
              #new_config = b"\x00\x00\x00\x00\x00\x01"
              dataId = b'\x01\x42'
              WRITE_DAT_REQUEST = b'\x2e'
              WRITE_DAT_RESPONSE = b'\x68'
              query = IsoTpParallelQuery(sendcan, logcan, CP.sccBus, [rdr_fw_address], [WRITE_DAT_REQUEST+dataId+new_config], [WRITE_DAT_RESPONSE], debug=True)
              query.get_data(0)
              print(f"Try {i+1}")
              break
            break
          except Exception as e:
            print(f"Failed {i}: {e}") 
      except Exception as e:
        print("Failed to enable tracks" + str(e))
  print("END Try to enable radar tracks")
  # END try to enable radar tracks

