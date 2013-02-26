from django.db import models
from django.db.models import Q

class ExpTagTable(models.Model):
    tagId = models.AutoField(primary_key=True)
    tagName = models.CharField(max_length=200)
    tagDesc = models.CharField(max_length=200)
    
    def __unicode__(self):
        return self.tagName
    
    class Meta:
        db_table = 'tagTable'


class Experiment(models.Model):
    idexperiment = models.AutoField(primary_key=True)
    experimentName = models.CharField(max_length=200)
    date = models.DateField()
    mask = models.BooleanField()
    tag = models.ManyToManyField(ExpTagTable,through='ExpToTag')
    
    def __unicode__(self):
        return self.experimentName
    
    class Meta:
        db_table = 'experiment'
   
   
class ExpToTag(models.Model):
    idexperimentTag = models.AutoField(primary_key=True)
    tagTable_tagId = models.ForeignKey(ExpTagTable,db_column='tagTable_tagId')
    experiment_idexperiment = models.ForeignKey(Experiment,db_column='experiment_idexperiment')
    
    class Meta:
        db_table = 'experimentTags'
 
 
class Record(models.Model):
    idrecord = models.AutoField(primary_key=True)
    experiment_idexperiment = models.ForeignKey(Experiment,db_column='experiment_idexperiment')
    filename = models.CharField(max_length=200)
    time = models.TimeField()
    numberofresponses = models.IntegerField()
    hardError = models.BooleanField()
    softError = models.BooleanField()
    
    def __unicode__(self):
        return self.filename
    
    class Meta:
        db_table = 'record'
   
   
class RecordTagTable(models.Model):
    idrecordTags = models.AutoField(primary_key=True)
    tagName = models.CharField(max_length=200)
    tagDesc = models.CharField(max_length=200)
    
    def __unicode__(self):
        return self.tagName
    
    class Meta:
        db_table = 'recordTags'
    
    
class RecordToTag(models.Model):
    idrecordToTags = models.AutoField(primary_key=True)
    recordTags_idrecordTags = models.ForeignKey(RecordTagTable,db_column='recordTags_idrecordTags')
    record_idrecord = models.ForeignKey(Record,db_column='record_idrecord')
    
    class Meta:
        db_table = 'recordToTags'
        

class SignalProp(models.Model):
    idsignalProperties = models.AutoField(primary_key=True)
    snr = models.FloatField()
    std = models.FloatField()
    mainLevel = models.IntegerField()
    ptp = models.FloatField()
    record_idrecord = models.ForeignKey(Record,db_column='record_idrecord')
    
    def __unicode__(self):
        return str(self.snr)
    
    class Meta:
        db_table = 'signalProperties'


class Response(models.Model):
    idresponses = models.AutoField(primary_key=True)    
    record_idrecord = models.ForeignKey(Record,db_column='record_idrecord')
    number = models.IntegerField()
    numberofspikes = models.IntegerField()
    length = models.FloatField()
    vpsp = models.FloatField()
    epspFront = models.FloatField()
    epspBack = models.FloatField()
    epileptStd = models.FloatField()
    epspArea = models.FloatField()
    
    def __unicode__(self):
        return str(self.number)
    
    class Meta:
        db_table = 'responses'
   
   
class Spike(models.Model):
    idspikes = models.AutoField(primary_key=True)    
    responses_idresponses = models.ForeignKey(Response,db_column='responses_idresponses')
    ampl = models.FloatField()
    number = models.IntegerField()
    length = models.FloatField()
    maxDiff = models.FloatField()
    angle1 = models.FloatField()
    angle2 = models.FloatField()
    delay = models.FloatField()
    maxtomin = models.FloatField()
    area = models.FloatField()
    fibre = models.BooleanField()
    manual = models.BooleanField()
            
    class Meta:
        db_table = 'spikes'
        
    def __unicode__(self):
        return str(self.idspikes)