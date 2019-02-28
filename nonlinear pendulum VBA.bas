Attribute VB_Name = "Module1"
Public Function corrector_predictor(qi, qdi, dt) As Double()

    Dim return_(1) As Double

    qs = qi + dt * qdi
    qds = qdi - dt * Sin(qi)
    
    qii = qi + dt / 2 * (qds + qdi)
    qdii = qdi - dt / 2 * (Sin(qs) + Sin(qi))
    
    return_(0) = qii
    return_(1) = qdii
    
    corrector_predictor = return_
    
End Function


Public Function rk_4(qi, qdi, dt) As Double()

    Dim return_(1) As Double
    
    '   dq/dt=F(qd)=qd
    '   d qd/dt =G(q)=-sin(q)
    

     qs = qi + dt / 2 * qdi
    qds = qdi - dt / 2 * Sin(qi)
    
     qss = qi + dt / 2 * qds
    qdss = qdi - dt / 2 * Sin(qs)
    
     qsss = qi + dt * qdss
    qdsss = qdi - dt * Sin(qss)
    
     qii = qi + dt / 6 * (qdi + 2 * qds + 2 * qdss + qdss)
    qdii = qdi - dt / 6 * (Sin(qi) + 2 * Sin(qs) + 2 * Sin(qss) + Sin(qsss))
    

    return_(0) = qii
    return_(1) = qdii
    
    rk_4 = return_
    
End Function


Public Function explicit(qi, qdi, dt) As Double()

    Dim return_(1) As Double

    qii = qi + dt * qdi
    qdii = qdi - dt * Sin(qi)
    
    return_(0) = qii
    return_(1) = qdii
    
    explicit = return_
    
End Function

Public Function midpoint(qi, qdi, dt) As Double()

    Dim return_(1) As Double

    qs = qi
    
    For i = 1 To 10
        qii = (Cos(qs) * dt ^ 2 * qs - Sin(qi) * dt ^ 2 - Sin(qs) * dt ^ 2 + 4# * qdi * dt + 4# * qi) / (dt ^ 2 * Cos(qs) + 4#)
        qdii = -(dt ^ 2 * Cos(qs) * qdi + 2# * dt * Cos(qs) * qi - 2# * Cos(qs) * dt * qs + 2# * Sin(qi) * dt + 2# * Sin(qs) * dt - 4# * qdi) / (dt ^ 2 * Cos(qs) + 4#)
        qs = qii
    Next
    
    return_(0) = qii
    return_(1) = qdii
    
    midpoint = return_
    
End Function

