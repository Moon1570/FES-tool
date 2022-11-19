from django.db import models

# Create your models here.
class Patient(models.Model):
    name = models.CharField(max_length=100)
    weight = models.IntegerField()
    intvalTime =  models.IntegerField(default=14)
    
    
    