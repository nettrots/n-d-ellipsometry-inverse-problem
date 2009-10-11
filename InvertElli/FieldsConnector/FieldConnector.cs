using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    //C# 3.0 inferring sugar
    public static class FieldConnector
    {
        public static FieldConnector<TObjectInterface, TObjectFied> Tie<TObjectInterface, TObjectFied>(TObjectInterface interface_field, TObjectFied internal_field, string interfacePr, string fieldPr, string name)
        {
            return new FieldConnector<TObjectInterface, TObjectFied>(interface_field, interfacePr, internal_field, fieldPr, name);
        }
 
    }
    public class FieldConnector< TObjectInterface,TObjectFied> : AbstractFieldConnector
    {
        private string _fieldPr;
        private string _interfacePr;

        public FieldConnector(TObjectInterface interface_field, string interfacePr, TObjectFied internal_field, string fieldPr, string name)
        {
            _field = internal_field;
            _interface = interface_field;
            _interfacePr = interfacePr;
            _fieldPr = fieldPr;
            Name = name;
        }
       
        public override void syncTo()
        {
            //TODO: automatic type conversation Tostring

            object obj = Convert.ToString(typeof (TObjectFied).GetProperty(_fieldPr).GetValue(_field, null));
            typeof(TObjectInterface).GetProperty(_interfacePr).SetValue(_interface, obj, null);
        }

        public override void syncFrom()
        {
            
            object obj = typeof(TObjectInterface).GetProperty(_interfacePr).GetValue(_interface, null);
            Type t = typeof (TObjectFied).GetProperty(_fieldPr).PropertyType;
            if(t==Type.GetType("System.Double"))
                obj = Convert.ToDouble(obj);
            else
            if (t == Type.GetType("System.Int32"))
                obj = Convert.ToInt32(obj);
            else
                throw new Exception("Can't convert interface to field");
            typeof(TObjectFied).GetProperty(_fieldPr).SetValue(_field, obj, null);
        }

    }
}
