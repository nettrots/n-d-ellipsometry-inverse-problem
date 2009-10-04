using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    public class FieldConnector<TObjectFied, TObjectInterface> : AbstractFieldConnector
    {
        private TObjectFied _field;
        private TObjectInterface _interface;
        private string _fieldPr;
        private string _interfacePr;

        FieldConnector(TObjectInterface interface_field, string interfacePr, TObjectFied internal_field, string fieldPr, string name)
        {
            _field = internal_field;
            _interface = interface_field;
            _interfacePr = interfacePr;
            _fieldPr = fieldPr;
            Name = name;
        }

        protected override void syncTo()
        {
            object obj = typeof(TObjectFied).GetProperty(_fieldPr).GetValue(_field, null);
            typeof(TObjectInterface).GetProperty(_interfacePr).SetValue(_interface, obj, null);
        }

        protected override void syncFrom()
        {
            object obj = typeof(TObjectInterface).GetProperty(_interfacePr).GetValue(_interface, null);
            typeof(TObjectFied).GetProperty(_fieldPr).SetValue(_field, obj, null);
        }

    }
}
